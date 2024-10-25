#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz, coo_matrix, diags, isspmatrix_csr, spdiags
import statsmodels.api as sm
from collections.abc import Iterable

# Define the standardize function outside the class
def standardize(array):
    std = np.std(array)
    if std == 0:
        return np.zeros_like(array)
    else:
        return (array - np.mean(array)) / std

class Normalization:
    def __init__(self):
        pass

    def preprocess(self, contig_file, contact_matrix_file, output_path, min_len=1000, min_signal=2, thres=5):
        self.min_len = min_len
        self.min_signal = min_signal
        self.output_path = output_path
        self.thres = thres

        # Load contig information
        names = ['contig_name', 'sites', 'length', 'coverage']
        self.contig_info = pd.read_csv(
            contig_file,
            header=None,
            names=names,
            dtype={'sites': float, 'length': float, 'coverage': float}
        )

        # Filter contigs based on min_len
        self.contig_info = self.contig_info[self.contig_info['length'] >= self.min_len].reset_index(drop=True)

        # Load contact matrix
        contact_matrix_full = load_npz(contact_matrix_file).tocoo()

        # No need to reindex if indices match
        self.contact_matrix = contact_matrix_full

        # Make output folder
        os.makedirs(self.output_path, exist_ok=True)


    def raw(self):
        """
        Do not conduct any normalization, just denoise.
        """
        try:
            contact_matrix_raw = self.contact_matrix.copy()
            self.denoise(contact_matrix_raw, 'raw')
            del contact_matrix_raw
        except Exception as e:
            print(f"Error during raw normalization: {e}")

    def normcc(self, epsilon=1):
        try:
            contact_matrix = self.contact_matrix.copy()
            covcc = contact_matrix.diagonal()
            contact_matrix.setdiag(0)
            signal = contact_matrix.max(axis=1).toarray().ravel()
            site = self.contig_info['sites'].values
            length = self.contig_info['length'].values

            # Debugging statements
            print(f"Length of site: {len(site)}")
            print(f"Length of length: {len(length)}")
            print(f"Length of covcc: {len(covcc)}")
            print(f"Length of signal: {len(signal)}")

            contig_info_normcc = pd.DataFrame({
                'site': site,
                'length': length,
                'covcc': covcc,
                'signal': signal
            })

            contig_info_normcc['sample_site'] = np.log(contig_info_normcc['site'] + epsilon)
            contig_info_normcc['sample_len'] = np.log(contig_info_normcc['length'])
            contig_info_normcc['sample_covcc'] = np.log(contig_info_normcc['covcc'] + epsilon)
            exog = contig_info_normcc[['sample_site', 'sample_len', 'sample_covcc']]
            endog = contig_info_normcc['signal']
            exog = sm.add_constant(exog)
            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit(method="lbfgs")
            norm_result = res.params.values

            if norm_result is None:
                print("normCC normalization failed.")
                return None

            linear_predictor = np.dot(exog, norm_result)
            expected_signal = np.exp(linear_predictor)
            scal = np.max(expected_signal)

            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(contact_matrix.row, contact_matrix.col, contact_matrix.data):
                mu_i = expected_signal[i]
                mu_j = expected_signal[j]
                normalized_value = scal * v / np.sqrt(mu_i * mu_j)
                normalized_data.append(normalized_value)

            # Create normalized contact matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (contact_matrix.row, contact_matrix.col)),
                shape=contact_matrix.shape
            )

            self.denoise(normalized_contact_matrix, 'normcc')
            del contact_matrix, normalized_contact_matrix

        except Exception as e:
            print(f"Error during normCC normalization: {e}")
            return None

    def hiczin(self, epsilon=1):
        try:
            contact_matrix = self.contact_matrix.copy()
            contact_matrix.setdiag(0)

            # Log transformation
            contig_info_hiczin = self.contig_info.copy()
            contig_info_hiczin['site'] = contig_info_hiczin['sites'] + epsilon

            # Handle missing or zero coverage
            if 'coverage' not in contig_info_hiczin.columns or contig_info_hiczin['coverage'].isnull().any():
                print("Warning: 'coverage' column is missing or contains NaNs. Replacing with epsilon.")
                contig_info_hiczin['coverage'] = epsilon
            else:
                min_non_zero = contig_info_hiczin['coverage'][contig_info_hiczin['coverage'] > 0].min()
                if pd.isna(min_non_zero):
                    min_non_zero = epsilon
                contig_info_hiczin['coverage'] = contig_info_hiczin['coverage'].replace(0, min_non_zero)

            map_x = contact_matrix.row
            map_y = contact_matrix.col
            map_data = contact_matrix.data
            index = map_x < map_y
            map_x = map_x[index]
            map_y = map_y[index]
            map_data = map_data[index]

            sample_site = np.log(contig_info_hiczin['site'].values[map_x] * contig_info_hiczin['site'].values[map_y])
            sample_len = np.log(contig_info_hiczin['length'].values[map_x] * contig_info_hiczin['length'].values[map_y])
            sample_cov = np.log(contig_info_hiczin['coverage'].values[map_x] * contig_info_hiczin['coverage'].values[map_y])

            # Use the safe standardization function
            sample_site = standardize(sample_site)
            sample_len = standardize(sample_len)
            sample_cov = standardize(sample_cov)

            # Debugging statements
            print(f"Standard deviation of sample_site: {np.std(sample_site)}")
            print(f"Standard deviation of sample_len: {np.std(sample_len)}")
            print(f"Standard deviation of sample_cov: {np.std(sample_cov)}")

            data_hiczin = pd.DataFrame({
                'sample_site': sample_site,
                'sample_len': sample_len,
                'sample_cov': sample_cov,
                'sampleCon': map_data
            })

            exog = data_hiczin[['sample_site', 'sample_len', 'sample_cov']]
            endog = data_hiczin['sampleCon']

            exog = sm.add_constant(exog)

            # Check for NaNs or Infs
            print("Checking exog for NaNs or Infs")
            print(f"exog contains NaNs: {np.isnan(exog).any()}")
            print(f"exog contains Infs: {np.isinf(exog).any()}")

            # Modified condition
            if np.isnan(exog.values).any() or np.isinf(exog.values).any():
                raise ValueError("exog contains inf or nans")

            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit()
            norm_result = res.params.values
            linear_predictor = np.dot(exog, norm_result)
            expected_signal = np.exp(linear_predictor)
            normalized_data = map_data / expected_signal

            # Create normalized contact matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (map_x, map_y)),
                shape=contact_matrix.shape
            )
            normalized_contact_matrix = normalized_contact_matrix + normalized_contact_matrix.transpose()

            self.denoise(normalized_contact_matrix, 'hiczin')
            del contact_matrix, normalized_contact_matrix

        except Exception as e:
            print(f"Error during HiCzin normalization: {e}")
            return None



    def bin3c(self, epsilon=1, max_iter=1000, tol=1e-6):
        try:
            # Get number of sites, avoid division by zero
            num_sites = self.contig_info['sites'].values
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
            bistochastic_matrix, _ = self._bisto_seq(normalized_contact_matrix, max_iter, tol)

            self.denoise(bistochastic_matrix, 'bin3c')

        except Exception as e:
            print(f"Error during bin3C normalization: {e}")
            return None


    def metator(self, epsilon=1):
        """
        Perform MetaTOR normalization.
        """
        try:
            covcc = self.contact_matrix.tocsr().diagonal() #MetaTOR use coverage instead of within-contig contacts (i.e., covcc)
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
        try:
            # Ensure the matrix is in COO format
            if not isinstance(_norm_matrix, coo_matrix):
                _norm_matrix = _norm_matrix.tocoo()

            denoised_matrix_file = os.path.join(self.output_path, f'denoised_contact_matrix_{suffix}.npz')

            if _norm_matrix.nnz == 0:
                print(f"Warning: The contact matrix '{suffix}' is empty. Skipping denoising.")
                # Save an empty matrix
                denoised_contact_matrix = coo_matrix(_norm_matrix.shape)
                save_npz(denoised_matrix_file, denoised_contact_matrix)
                print(f"Empty denoised contact matrix saved to {denoised_matrix_file}")
                return denoised_contact_matrix

            # Calculate threshold to discard lowest p percent
            if self.thres <= 0 or self.thres >= 100:
                print("Error: The threshold percentage for spurious contact detection must be between 0 and 100.")
                return None

            threshold = np.percentile(_norm_matrix.data, self.thres)
            mask = _norm_matrix.data > threshold

            if not np.any(mask):
                print(f"Warning: No contacts exceed the threshold in '{suffix}'. Denoised matrix will be empty.")
                denoised_contact_matrix = coo_matrix(_norm_matrix.shape)
            else:
                denoised_contact_matrix = coo_matrix(
                    (_norm_matrix.data[mask], (_norm_matrix.row[mask], _norm_matrix.col[mask])),
                    shape=_norm_matrix.shape
                )

            # Save denoised contact matrix
            save_npz(denoised_matrix_file, denoised_contact_matrix)
            print(f"Denoised normalized contact matrix '{suffix}' saved to {denoised_matrix_file}")

            return denoised_contact_matrix
        except Exception as e:
            print(f"Error during denoising normalization for '{suffix}': {e}")
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
        if x0 is None:
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
                Z = rk / v
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
    parser = argparse.ArgumentParser(description="Normalization tool for MetaHit pipeline.")
    subparsers = parser.add_subparsers(dest='command', help='Available normalization commands')

    # Raw normalization
    parser_raw = subparsers.add_parser('raw', help='Perform raw normalization')
    parser_raw.add_argument('--contig_file', '-c', required=True, help='Path to contig_info.csv')
    parser_raw.add_argument('--contact_matrix_file', '-m', required=True, help='Path to contact_matrix.npz')
    parser_raw.add_argument('--output_path', '-o', required=True, help='Output directory')
    parser_raw.add_argument('--min_len', type=int, default=1000, help='Minimum contig length')
    parser_raw.add_argument('--min_signal', type=int, default=2, help='Minimum signal')
    parser_raw.add_argument('--thres', type=float, default=5, help='Threshold percentage for denoising (0-100)')

    # normCC normalization
    parser_normcc = subparsers.add_parser('normcc', help='Perform normCC normalization')
    parser_normcc.add_argument('--contig_file', '-c', required=True, help='Path to contig_info.csv')
    parser_normcc.add_argument('--contact_matrix_file', '-m', required=True, help='Path to contact_matrix.npz')
    parser_normcc.add_argument('--output_path', '-o', required=True, help='Output directory')
    parser_normcc.add_argument('--thres', type=float, default=5, help='Threshold percentage for denoising (0-100)')


    # Subparser for HiCzin normalization
    parser_hiczin = subparsers.add_parser('hiczin', help='Perform HiCzin normalization')
    parser_hiczin.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_hiczin.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_hiczin.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_hiczin.add_argument('--epsilon', type=float, default=1, help='Epsilon value')
    parser_hiczin.add_argument('--thres', type=float, default=5, help='Threshold percentage for denoising (0-100)')

    # Subparser for bin3C normalization
    parser_bin3c = subparsers.add_parser('bin3c', help='Perform bin3C normalization')
    parser_bin3c.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_bin3c.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_bin3c.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_bin3c.add_argument('--epsilon', type=float, default=1, help='Epsilon value')
    parser_bin3c.add_argument('--max_iter', type=int, default=1000, help='Maximum iterations for Sinkhorn-Knopp')
    parser_bin3c.add_argument('--tol', type=float, default=1e-6, help='Tolerance for convergence')
    parser_bin3c.add_argument('--thres', type=float, default=5, help='Threshold percentage for denoising (0-100)')

    # Subparser for MetaTOR normalization
    parser_metator = subparsers.add_parser('metator', help='Perform MetaTOR normalization')
    parser_metator.add_argument('-c', '--contig_file', required=True, help='Path to contig_info.csv')
    parser_metator.add_argument('-m', '--contact_matrix_file', required=True, help='Path to contact_matrix.npz')
    parser_metator.add_argument('-o', '--output_path', required=True, help='Output directory')
    parser_metator.add_argument('--epsilon', type=float, default=1, help='Epsilon value')
    parser_metator.add_argument('--thres', type=float, default=5, help='Threshold percentage for denoising (0-100)')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    normalizer = Normalization()

    if args.command == 'raw':
        normalizer.preprocess(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            min_len=args.min_len,
            min_signal=args.min_signal,
            thres=args.thres
        )
        normalizer.raw()
    elif args.command == 'normcc':
        normalizer.preprocess(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            thres=args.thres
        )
        normalizer.normcc()

    elif args.command == 'hiczin':
        normalizer.preprocess(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            thres=args.thres
        )
        normalizer.hiczin(epsilon=args.epsilon)

    elif args.command == 'bin3c':
        normalizer.preprocess(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            thres=args.thres
        )
        normalizer.bin3c(epsilon=args.epsilon, max_iter=args.max_iter, tol=args.tol)

    elif args.command == 'metator':
        normalizer.preprocess(
            contig_file=args.contig_file,
            contact_matrix_file=args.contact_matrix_file,
            output_path=args.output_path,
            thres=args.thres
        )
        normalizer.metator(epsilon=args.epsilon)

if __name__ == "__main__":
    main()
