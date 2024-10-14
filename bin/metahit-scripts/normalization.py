#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz, coo_matrix, diags, isspmatrix_csr, spdiags
import statsmodels.api as sm
from collections import namedtuple
from collections.abc import Iterable

# Ensure the script directory is in sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))


"""
Input:
1) Raw Hi-C contact matrix: e.g. Raw_contact_matrix.npz
2) Contig info file: e.g. csv four columns, no header name: contig name; number of sites; contig length; contig coverage
3) parameter: threshold of spurious contacts
By default, we will do spurious contact detection within the normalization steps
contig info can be three columns or four columns depending on whether the depth is computed or not
if four columns then depth is computed
"""

class Normalization:
    def __init__(self):
        pass

    def preprocess(self, contig_file, raw_contact_file, output_path, min_len=1000, min_signal=2, thres = 5):
        self.min_len = min_len
        self.min_signal = min_signal
        self.contact_matrix = load_npz(raw_contact_file).tocoo()
        names = ['contig_name', 'site', 'length', 'coverage']
        self.contig_info = pd.read_csv(
                contig_file,
                header=None,
                names=names,
                dtype={'site': float, 'length': float, 'coverage': float}
            )
        #####Do contig filtering according to min len and min signal
        _m = self.contact_matrix
        _m = _m.tolil(True)
        _diag = _m.tocsr().diagonal()
        _m.setdiag(0)
        _sig = np.asarray(_m.tocsr().max(axis=0).todense()).ravel()
        _contig_id = []
        for i in range(_m.shape[0]):
            if _sig[i] >= self.min_signal and self.contig_info['sites'].iloc[i]>0:
                _contig_id.append(i)
        del _m, _diag, _sig
        
        ######filter contig info and raw Hi-C matrix
        self.contig_info = self.contig_info.iloc[_contig_id]

        self.contact_matrix = self.contact_matrix.tocsr()
        self.contact_matrix = self.contact_matrix[_contig_id , :]
        self.contact_matrix = self.contact_matrix.tocsc()
        self.contact_matrix = self.contact_matrix[: , _contig_id]
        self.contact_matrix = self.contact_matrix.tocoo()
        del contig_id
        
        self.output_path = output_path
        self.thres = thres
        
        # Make output folder 
        os.makedirs(self.output_path, exist_ok=True)

    def raw(self, min_len=1000, min_signal=2):
        """
        Do not conduct any normalization
        """
        try:
            contact_matrix_raw = self.contact_matrix.copy()
            self.denoise(contact_matrix_raw, 'raw')
            del contact_matrix_raw

        except Exception as e:
            print(f"Error during raw normalization: {e}")
            return None

    def normcc(self, epsilon=1):
        """
        Perform normcc normalization using the custom normcc_local function.
        """
        try:
            contact_matrix = self.contact_matrix.copy()
            covcc = contact_matrix.tocsr().diagonal()
            contact_matrix.setdiag(0)
            signal = np.asarray(contact_matrix.tocsr().max(axis=0).todense()).ravel()
            site = self.contig_info['site'].values  
            length = self.contig_info['length'].values 

            contig_info_normcc = pd.DataFrame({
                                                'site': site,
                                                'length': length,
                                                'covcc': covcc,
                                                'signal': signal
                                            })
            del site, length, covcc, signal
            
            contig_info_normcc['sample_site'] = np.log(contig_info_normcc['site'] + epsilon)
            contig_info_normcc['sample_len'] = np.log(contig_info_normcc['length'])
            contig_info_normcc['sample_covcc'] = np.log(contig_info_normcc['covcc']+ epsilon)
            exog = contig_info_normcc[['sample_site', 'sample_len', 'sample_covcc']]
            endog = contig_info_normcc[['signal']]
            exog = sm.add_constant(exog)
            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit(method="lbfgs")
            norm_result = res.params.to_list()
            

            if norm_result is None:
                print("normcc normalization failed.")
                return None

            linear_predictor = np.dot(exog, norm_result)
            expected_signal = np.exp(linear_predictor)
            scal = np.max(expected_signal)

            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(contact_matrix.row, contact_matrix.col, contact_matrix.data):
                mu_i = expected_signal[i]
                mu_j = expected_signal[j]
                normalized_value = scal*v/np.sqrt(mu_i*mu_j)
                normalized_data.append(normalized_value)

            # Create normalized contact matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (contact_matrix.row, contact_matrix.col)),
                shape=contact_matrix.shape
            )
        
            self.denoise(normalized_contact_matrix, 'normcc')
            del contact_matrix, normalized_contact_matrix

        except Exception as e:
            print(f"Error during normcc normalization: {e}")
            return None
        

    def hiczin(self, epsilon=1):
        """
        hiczin require coverage information, if there is no coverage information, report error
        """
        try:
            contact_matrix = self.contact_matrix.copy()
            contact_matrix.setdiag(0)
            # Log transformation
            contig_info_hiczin = self.contig_info
            contig_info_hiczin['site'] = contig_info_hiczin['site'] + epsilon
            #Replace the zero coverage (which is very rare) to the min non zero entry
            min_non_zero = contig_info_hiczin['coverage'][contig_info_hiczin['coverage'] > 0].min()
            contig_info_hiczin['coverage'] = contig_info_hiczin['coverage'].replace(0, min_non_zero)    
            map_x = self.contact_matrix.row
            map_y = self.contact_matrix.col
            map_data = self.contact_matrix.data
            index = map_x < map_y
            map_x = map_x[index]
            map_y = map_y[index]
            map_data = map_data[index]
            sample_len = np.zeros(len(map_x))
            sample_site = np.zeros(len(map_x))
            sample_cov = np.zeros(len(map_x))           
            for i in range(len(map_x)):
                idx1 = int(map_x[i])
                idx2 = int(map_y[i])

                sample_site[i] = np.log(contig_info_hiczin.iloc[idx1]['site'] * contig_info_hiczin.iloc[idx2]['site'])
                sample_len[i] = np.log(contig_info_hiczin.iloc[idx1]['length'] * contig_info_hiczin.iloc[idx2]['length'])
                sample_cov[i] = np.log(contig_info_hiczin.iloc[idx1]['coverage'] * contig_info_hiczin.iloc[idx2]['coverage'])
            
            sample_site = (sample_site - np.mean(sample_site)) / np.std(sample_site)
            sample_len = (sample_len - np.mean(sample_len)) / np.std(sample_len)
            sample_cov = (sample_cov - np.mean(sample_cov)) / np.std(sample_cov)     
            data_hiczin = pd.DataFrame({
                                        'sample_site': sample_site,
                                        'sample_len': sample_len,
                                        'sample_cov': sample_cov,
                                        'sampleCon': map_data  
                                        })
            

            exog = data_hiczin[['sample_site', 'sample_len', 'sample_cov']]
            endog = data_hiczin['sampleCon']

            exog = sm.add_constant(exog)
            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit()
            norm_result = res.params.to_list()
            linear_predictor = np.dot(exog, norm_result)
            expected_signal = np.exp(linear_predictor)
            normalized_data = map_data/expected_signal
            # Create normalized contact matrix, it initilize as a lower triangle matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (map_x, map_y)),
                shape=contact_matrix.shape
            )
            normalized_contact_matrix = normalized_contact_matrix + normalized_contact_matrix.transpose()
        
            self.denoise(normalized_contact_matrix, 'hiczin')  
            del contact_matrix, normalized_contact_matrix       

        except Exception as e:
            print(f"Error during HiCzin parameter calculation: {e}")
            return None


    def bin3c(self, epsilon=1, max_iter=1000, tol=1e-6):
        """
        Perform bin3C normalization.
        """
        try:
            # Get number of sites, avoid division by zero
            num_sites = self.contig_info['site'].values
            num_sites = num_sites + epsilon

            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(self.contact_matrix.row, self.contact_matrix.col, self.contact_matrix.data):
                s_i = num_sites[i]
                s_j = num_sites[j]
                norm_value = v / (s_i * s_j)
                normalized_data.append(norm_value)

            normalized_contact_matrix = coo_matrix(
                (normalized_data, (self.contact_matrix.row, self.contact_matrix.col)),
                shape=self.contact_matrix.shape
            )

            # Apply KR algorithm
            bistochastic_matrix = self._bisto_seq(normalized_contact_matrix, max_iter, tol)

            self.denoise(bistochastic_matrix)
            
        except Exception as e:
            print(f"Error during bin3C normalization: {e}")
            return None


    def metator(self, epsilon=1):
        """
        Perform MetaTOR normalization.
        """
        try:
            covcc = self.contact_matrix.tocsr().diagonal()
            covcc = covcc + epsilon
            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(self.contact_matrix.row, self.contact_matrix.col, self.contact_matrix.data):
                cov_i = covcc[i]
                cov_j = covcc[j]
                norm_factor = np.sqrt(cov_i * cov_j)
                normalized_value = v / norm_factor
                normalized_data.append(normalized_value)

            # Create normalized contact matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (self.contact_matrix.row, self.contact_matrix.col)),
                shape=self.contact_matrix.shape
            )

            self.denoise(normalized_contact_matrix, 'metator')
            
        except Exception as e:
            print(f"Error during MetaTOR normalization: {e}")
            return None
        
    def denoise(self, _norm_matrix, suffix):
        """
        Remove spurious Hi-C contacts by discarding the lowest p percent of normalized contacts.
        Saves the denoised matrix as 'denoised_normalized_matrix.npz'.
        
        Parameters:
        - contact_matrix_file: Path to the normalized contact matrix (.npz).
        - output_path: Directory to save the denoised matrix.
        - p: Fraction of contacts to discard (default is 0.1 for 10%).
        """
        try:
            denoised_matrix_file = os.path.join(self.output_path, 'denoised_contact_matrix_' + str(suffix) + '.npz')

            if _norm_matrix == 0:
                print("Warning: The contact matrix is empty. Skipping denoising.")
                # Save an empty matrix
                denoised_contact_matrix = coo_matrix(_norm_matrix.shape)
                save_npz(denoised_matrix_file, denoised_contact_matrix)
                print(f"Empty denoised contact matrix saved to {denoised_matrix_file}")
                return denoised_contact_matrix


            # Calculate threshold to discard lowest p percent
            if self.thres <= 0 or self.thres >= 1:
                print("Error: The threshold percentage for spurious contact detection must be between 0 and 1.")
                return None

            threshold = np.percentile(_norm_matrix.data, self.thres * 100)
            mask = _norm_matrix.data > threshold

            if not np.any(mask):
                print("Warning: No contacts exceed the threshold. Denoised matrix will be empty.")
                denoised_contact_matrix = coo_matrix(_norm_matrix.shape)
            else:
                denoised_contact_matrix = coo_matrix(
                    (_norm_matrix.data[mask], (_norm_matrix.row[mask], _norm_matrix.col[mask])),
                    shape=_norm_matrix.shape
                )

            # Save denoised contact matrix
            save_npz(denoised_matrix_file, denoised_contact_matrix)
            print(f"Denoised normalized contact matrix saved to {denoised_matrix_file}")

            return denoised_contact_matrix
        except Exception as e:
            print(f"Error during denoising normalization: {e}")
            return None
        

    def _bisto_seq(self, m , max_iter, tol, x0=None, delta=0.1, Delta=3):
        _orig = m.copy()
        # replace 0 diagonals with 1, on the working matrix. This avoids potentially
        # exploding scale-factors. KR should be regularlized!
        m = m.tolil()
        is_zero = m.diagonal() == 0
        if np.any(is_zero):
            print('treating {} zeros on diagonal as ones'.format(is_zero.sum()))
            ix = np.where(is_zero)
            m[ix, ix] = 1
        if not isspmatrix_csr(m):
            m = m.tocsr()
        n = m.shape[0]
        e = np.ones(n)
        if not x0:
            x0 = e.copy()
        g = 0.9
        etamax = 0.1
        eta = etamax
        stop_tol = tol * 0.5
        x = x0.copy()
        rt = tol ** 2
        v = x * m.dot(x)
        rk = 1 - v
        rho_km1 = rk.T.dot(rk)  # transpose possibly implicit
        rout = rho_km1
        rold = rout
        n_iter = 0
        i = 0
        y = np.empty_like(e)
        while rout > rt and n_iter < max_iter:
            i += 1
            k = 0
            y[:] = e
            inner_tol = max(rout * eta ** 2, rt)
            while rho_km1 > inner_tol:
                k += 1
                if k == 1:
                    Z = rk / v
                    p = Z
                    rho_km1 = rk.T.dot(Z)
                else:
                    beta = rho_km1 / rho_km2
                    p = Z + beta * p
                w = x * m.dot(x * p) + v * p
                alpha = rho_km1 / p.T.dot(w)
                ap = alpha * p
                ynew = y + ap
                if np.amin(ynew) <= delta:
                    if delta == 0:
                        break
                    ind = np.where(ap < 0)[0]
                    gamma = np.amin((delta - y[ind]) / ap[ind])
                    y += gamma * ap
                    break
                if np.amax(ynew) >= Delta:
                    ind = np.where(ynew > Delta)[0]
                    gamma = np.amin((Delta - y[ind]) / ap[ind])
                    y += gamma * ap
                    break
                y = ynew
                rk = rk - alpha * w
                rho_km2 = rho_km1
                Z = rk * v
                rho_km1 = np.dot(rk.T, Z)
                if np.any(np.isnan(x)):
                    raise RuntimeError('scale vector has developed invalid values (NANs)!')
            x *= y
            v = x * m.dot(x)
            rk = 1 - v
            rho_km1 = np.dot(rk.T, rk)
            rout = rho_km1
            n_iter += k + 1
            rat = rout / rold
            rold = rout
            res_norm = np.sqrt(rout)
            eta_o = eta
            eta = g * rat
            if g * eta_o ** 2 > 0.1:
                eta = max(eta, g * eta_o ** 2)
            eta = max(min(eta, etamax), stop_tol / res_norm)
        if n_iter > max_iter:
            raise RuntimeError('matrix balancing failed to converge in {} iterations'.format(n_iter))
        del m
        print('It took {} iterations to achieve bistochasticity'.format(n_iter))
        if n_iter >= max_iter:
            print('Warning: maximum number of iterations ({}) reached without convergence'.format(max_iter))
        X = spdiags(x, 0, n, n, 'csr')
        return X.T.dot(_orig.dot(X)), x

def main():
    parser = argparse.ArgumentParser(
        description="Normalization tool for MetaHit pipeline."
    )
    subparsers = parser.add_subparsers(dest='command', help='Available normalization commands')

    # Subparser for raw normalization
    parser_raw = subparsers.add_parser('raw', help='Perform raw normalization')
    parser_raw.add_argument('-b', '--bam_file', required=True, help='Path to BAM file')
    parser_raw.add_argument('-f', '--fasta_file', required=True, help='Path to FASTA file')
    parser_raw.add_argument('-e', '--enzymes', nargs='+', required=True, help='List of enzymes')
    parser_raw.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_raw.add_argument('--min_mapq', type=int, default=30, help='Minimum MAPQ')
    parser_raw.add_argument('--min_len', type=int, default=1000, help='Minimum length')
    parser_raw.add_argument('--min_match', type=int, default=30, help='Minimum match')
    parser_raw.add_argument('--min_signal', type=int, default=2, help='Minimum signal')

    # Subparser for normcc_contact normalization
    parser_normcc = subparsers.add_parser('normcc_contact', help='Perform normcc_contact normalization')
    parser_normcc.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_normcc.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_normcc.add_argument('-o', '--output_path', required=True, help='Output directory')

    # Subparser for HiCzin normalization
    parser_hiczin = subparsers.add_parser('hiczin', help='Perform HiCzin normalization')
    parser_hiczin.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_hiczin.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_hiczin.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_hiczin.add_argument('--min_signal', type=int, default=2, help='Minimum signal')

    # Subparser for bin3c normalization
    parser_bin3c = subparsers.add_parser('bin3c', help='Perform bin3C normalization')
    parser_bin3c.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_bin3c.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_bin3c.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_bin3c.add_argument('--max_iter', type=int, default=1000, help='Maximum iterations for Sinkhorn-Knopp')
    parser_bin3c.add_argument('--tol', type=float, default=1e-6, help='Tolerance for convergence')

    # Subparser for MetaTOR normalization
    parser_metator = subparsers.add_parser('metator', help='Perform MetaTOR normalization')
    parser_metator.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_metator.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_metator.add_argument('-o', '--output_path', required=True, help='Output directory')

    # Subparser for denoising
    parser_denoise = subparsers.add_parser('denoise', help='Perform denoising on a contact matrix')
    parser_denoise.add_argument('-m', '--contact_matrix_file', required=True, help='Path to normalized_contact_matrix.npz')
    parser_denoise.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_denoise.add_argument('-p', '--p', type=float, default=0.1, help='Fraction of contacts to discard')
    parser_denoise.add_argument('--discard_file', help='Path to discard contacts CSV file')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    normalizer = Normalization()

    normalizer.preprocess(
        
    )

    if args.command == 'raw':
        cm = normalizer.raw(
            bam_file=args.bam_file,
            fasta_file=args.fasta_file,
            enzymes=args.enzymes,
            output_path=args.output_path,
            min_mapq=args.min_mapq,
            min_len=args.min_len,
            min_match=args.min_match,
            min_signal=args.min_signal
        )
        if cm is not None:
            print("Raw normalization completed successfully.")
        else:
            print("Raw normalization failed.")

    elif args.command == 'normcc_contact':
        norm_result = normalizer.normcc_contact(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path
        )
        if norm_result is not None:
            print("normCC normalization completed successfully.")
        else:
            print("normCC normalization failed.")

    elif args.command == 'hiczin':
        hiczin_map = normalizer.hiczin(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            min_signal=args.min_signal
        )
        if hiczin_map is not None:
            print("HiCzin normalization completed successfully.")
        else:
            print("HiCzin normalization failed.")

    elif args.command == 'bin3c':
        bistochastic_matrix = normalizer.bin3c(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            max_iter=args.max_iter,
            tol=args.tol
        )
        if bistochastic_matrix is not None:
            print("bin3C normalization completed successfully.")
        else:
            print("bin3C normalization failed.")

    elif args.command == 'metator':
        metator_matrix = normalizer.metator(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path
        )
        if metator_matrix is not None:
            print("MetaTOR normalization completed successfully.")
        else:
            print("MetaTOR normalization failed.")

    elif args.command == 'denoise':
        denoised_matrix = normalizer.denoise(
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            p=args.p,
            discard_file=args.discard_file
        )
        if denoised_matrix is not None:
            print("Denoising completed successfully.")
        else:
            print("Denoising failed.")

if __name__ == "__main__":
    main()
