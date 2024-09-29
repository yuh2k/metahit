import os
import sys
import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz, coo_matrix, diags
import statsmodels.api as sm
from collections import namedtuple

# Ensure the script directory is in sys.path
sys.path.append('bin/metahit-scripts')

from raw_contact import ContactMatrix as RawContactMatrix
from normcc import normcc_local
from hiczin_contact import HiCzinMap

class Normalization:
    def __init__(self):
        pass

    def raw(self, bam_file, fasta_file, enzymes, output_path,
            min_mapq=30, min_len=1000, min_match=30, min_signal=2):
        """
        Perform raw normalization using the ContactMatrix class from raw_contact.py.
        """
        print(f"Using raw normalization on BAM file: {bam_file} and FASTA file: {fasta_file}")
        try:
            os.makedirs(output_path, exist_ok=True)

            # Initialize ContactMatrix instance
            cm = RawContactMatrix(
                bam_file=bam_file,
                enzymes=enzymes,
                seq_file=fasta_file,
                path=output_path,
                min_mapq=min_mapq,
                min_len=min_len,
                min_match=min_match,
                min_signal=min_signal
            )

            # Save contact matrix
            matrix_file = os.path.join(output_path, 'contact_matrix.npz')
            save_npz(matrix_file, cm.seq_map)
            print(f"Contact matrix saved to {matrix_file}")

            # Save contig information
            contig_info_file = os.path.join(output_path, 'contig_info.csv')
            with open(contig_info_file, 'w') as f:
                for i, seq in enumerate(cm.seq_info):
                    f.write(f"{seq.name},{seq.sites},{seq.length},{seq.covcc},{cm.row_sum[i]}\n")
            print(f"Contig info saved to {contig_info_file}")

            print(f"Raw normalization result saved to {output_path}/contig_info.csv")
            return cm
        except Exception as e:
            print(f"Error during raw normalization: {e}")
            return None

    def normcc_contact(self, contig_file, contact_matrix_file, output_path):
        """
        Perform normcc normalization using the custom normcc_local function.
        """
        try:
            os.makedirs(output_path, exist_ok=True)

            coefficients_file = os.path.join(output_path, 'normcc_coefficients.csv')

            # Execute normcc_local to obtain normalization coefficients
            norm_result = normcc_local(contig_file, output_file=coefficients_file)

            if norm_result is None:
                print("normcc normalization failed.")
                return None

            print(f"normcc normalization coefficients saved to {coefficients_file}")
            print(f"Length of norm_result in normcc_contact: {len(norm_result)}")
            print(f"norm_result in normcc_contact: {norm_result}")

            # Load original contact matrix
            contact_matrix = load_npz(contact_matrix_file).tocoo()

            # Load contig information
            names = ['contig_name', 'site', 'length', 'covcc', 'signal']
            df = pd.read_csv(
                contig_file,
                header=None,
                names=names,
                dtype={'site': float, 'length': float, 'covcc': float, 'signal': float}
            )

            print("Data types of contig info:")
            print(df.dtypes)

            # Handle missing values
            if df[['site', 'length', 'covcc', 'signal']].isnull().values.any():
                print("Warning: Missing values found in 'site', 'length', 'covcc', or 'signal' columns.")
                df = df.dropna(subset=['site', 'length', 'covcc', 'signal'])

            # Calculate expected signal
            epsilon = 1e-6
            df['sample_site'] = np.log(df['site'] + epsilon)
            df['sample_len'] = np.log(df['length'] + epsilon)
            df['sample_covcc'] = np.log(df['covcc'] + epsilon)

            exog = df[['sample_site', 'sample_len', 'sample_covcc']]
            exog = sm.add_constant(exog)

            linear_predictor = np.dot(exog, norm_result)
            expected_signal = np.exp(linear_predictor)
            expected_signal += epsilon  # Prevent division by zero

            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(contact_matrix.row, contact_matrix.col, contact_matrix.data):
                mu_i = expected_signal[i]
                mu_j = expected_signal[j]
                norm_factor = (mu_i * mu_j) ** 0.5  # Geometric mean
                normalized_value = v / norm_factor
                normalized_data.append(normalized_value)

            # Create normalized contact matrix
            normalized_contact_matrix = coo_matrix(
                (normalized_data, (contact_matrix.row, contact_matrix.col)),
                shape=contact_matrix.shape
            )

            # Save normalized contact matrix
            normalized_matrix_file = os.path.join(output_path, 'normalized_contact_matrix.npz')
            save_npz(normalized_matrix_file, normalized_contact_matrix)
            print(f"Normalized contact matrix saved to {normalized_matrix_file}")

            return norm_result
        except Exception as e:
            print(f"Error during normcc normalization: {e}")
            return None

    def hiczin_normcc(self, contig_file, output_file=None):
        """
        Calculate normalization parameters required for HiCzin.
        """
        print(f"Calculating HiCzin normalization parameters using contig file: {contig_file}")
        try:
            names = ['contig_name', 'site', 'length', 'covcc', 'signal']
            df = pd.read_csv(
                contig_file,
                header=None,
                names=names,
                dtype={'site': float, 'length': float, 'covcc': float, 'signal': float}
            )

            if df[['site', 'length', 'covcc', 'signal']].isnull().values.any():
                print("Warning: Missing values found.")
                df = df.dropna(subset=['site', 'length', 'covcc', 'signal'])

            if df.shape[0] <= 4:
                print("Not enough data points to fit the model.")
                return None

            # Log transformation
            epsilon = 1e-6
            df['sample_site'] = np.log(df['site'] + epsilon)
            df['sample_len'] = np.log(df['length'] + epsilon)
            df['sample_covcc'] = np.log(df['covcc'] + epsilon)

            exog = df[['sample_site', 'sample_len', 'sample_covcc']]
            endog = df['signal']

            exog = sm.add_constant(exog)

            # Fit Negative Binomial model
            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial())
            res = glm_nb.fit()

            alpha = res.scale
            print(f"Estimated alpha for HiCzin: {alpha}")

            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=alpha))
            res = glm_nb.fit()

            norm_result = res.params.to_list()

            # Calculate threshold and statistics
            threshold = df['signal'][df['signal'] > 0].min()
            if pd.isna(threshold):
                threshold = epsilon

            s_mean = df['sample_site'].mean()
            s_std = df['sample_site'].std(ddof=0)

            l_mean = df['sample_len'].mean()
            l_std = df['sample_len'].std(ddof=0)

            c_mean = df['sample_covcc'].mean()
            c_std = df['sample_covcc'].std(ddof=0)

            norm_result.extend([threshold, s_mean, s_std, l_mean, l_std, c_mean, c_std])

            if output_file:
                output_dir = os.path.dirname(output_file)
                os.makedirs(output_dir, exist_ok=True)
                with open(output_file, 'w') as f:
                    f.write(','.join(map(str, norm_result)))
                    f.write('\n')
                print(f"HiCzin normalization parameters saved to {output_file}")

            return norm_result
        except Exception as e:
            print(f"Error during HiCzin parameter calculation: {e}")
            return None

    def hiczin(self, contig_file, contact_matrix_file, output_path, min_signal=2):
        """
        Perform HiCzin normalization.
        """
        try:
            os.makedirs(output_path, exist_ok=True)

            # Calculate normalization parameters for HiCzin
            hiczin_norm_result = self.hiczin_normcc(contig_file)
            if hiczin_norm_result is None:
                print("HiCzin normalization failed during parameter calculation.")
                return None

            # Load contig information
            names = ['contig_name', 'sites', 'length', 'covcc', 'signal']
            contig_info_df = pd.read_csv(
                contig_file,
                header=None,
                names=names,
                dtype={'sites': float, 'length': float, 'covcc': float, 'signal': float}
            )

            # Create contig_info list
            ContigInfo = namedtuple('ContigInfo', ['name', 'sites', 'length', 'cov', 'tax'])
            contig_info_list = []
            for idx, row in contig_info_df.iterrows():
                contig_info_list.append(
                    ContigInfo(
                        name=row['contig_name'],
                        sites=row['sites'],
                        length=row['length'],
                        cov=row['covcc'],
                        tax=''  # Leave empty if no tax information
                    )
                )

            # Load contact matrix
            seq_map = load_npz(contact_matrix_file)

            # Initialize HiCzinMap
            hiczin_map = HiCzinMap(
                path=output_path,
                contig_info=contig_info_list,
                seq_map=seq_map,
                norm_result=hiczin_norm_result,
                min_signal=min_signal
            )

            # Save normalized contact matrix
            hiczin_matrix_file = os.path.join(output_path, 'hiczin_contact_matrix.npz')
            hiczin_seq_map = hiczin_map.seq_map.tocoo()  # Convert to COO for saving
            save_npz(hiczin_matrix_file, hiczin_seq_map)
            print(f"HiCzin normalized contact matrix saved to {hiczin_matrix_file}")

            return hiczin_map
        except Exception as e:
            print(f"Error during HiCzin normalization: {e}")
            return None

    def bin3c(self, contig_file, contact_matrix_file, output_path, max_iter=1000, tol=1e-6):
        """
        Perform bin3C normalization.
        """
        try:
            os.makedirs(output_path, exist_ok=True)
            epsilon = 1e-10

            # Load contig information
            names = ['contig_name', 'site', 'length', 'covcc', 'signal']
            df = pd.read_csv(
                contig_file,
                header=None,
                names=names,
                dtype={'site': float, 'length': float, 'covcc': float, 'signal': float}
            )

            # Load raw contact matrix
            contact_matrix = load_npz(contact_matrix_file).tocoo()

            # Get number of sites, avoid division by zero
            num_sites = df['site'].values
            num_sites[num_sites == 0] = epsilon

            # Normalize contact values
            normalized_data = []
            for i, j, v in zip(contact_matrix.row, contact_matrix.col, contact_matrix.data):
                s_i = num_sites[i]
                s_j = num_sites[j]
                norm_value = v / (s_i * s_j)
                normalized_data.append(norm_value)

            normalized_contact_matrix = coo_matrix(
                (normalized_data, (contact_matrix.row, contact_matrix.col)),
                shape=contact_matrix.shape
            )

            # Apply Sinkhorn-Knopp algorithm
            bistochastic_matrix = self._bisto_seq(normalized_contact_matrix, max_iter, tol)

            # Save bistochastic matrix
            bin3c_matrix_file = os.path.join(output_path, 'bin3c_contact_matrix.npz')
            save_npz(bin3c_matrix_file, bistochastic_matrix)
            print(f"bin3C normalized contact matrix saved to {bin3c_matrix_file}")

            return bistochastic_matrix
        except Exception as e:
            print(f"Error during bin3C normalization: {e}")
            return None

    def _bisto_seq(self, matrix, max_iter=1000, tol=1e-6):
        """
        Convert matrix to bistochastic using Sinkhorn-Knopp algorithm.
        """
        epsilon = 1e-10
        A = matrix.tocsr()
        n = A.shape[0]
        r = np.ones(n)
        c = np.ones(n)

        for iter_num in range(max_iter):
            # Update row scaling factors
            row_sums = A.dot(c)
            row_sums[row_sums == 0] = epsilon
            r_new = 1.0 / row_sums
            r_change = np.linalg.norm(r_new - r)
            r = r_new

            # Update column scaling factors
            col_sums = A.T.dot(r)
            col_sums[col_sums == 0] = epsilon
            c_new = 1.0 / col_sums
            c_change = np.linalg.norm(c_new - c)
            c = c_new

            # Check for convergence
            if r_change < tol and c_change < tol:
                print(f"Sinkhorn-Knopp algorithm converged in {iter_num+1} iterations.")
                break
        else:
            print("Sinkhorn-Knopp algorithm did not converge within the maximum number of iterations.")

        # Create diagonal matrices for scaling
        D_r = diags(r)
        D_c = diags(c)

        # Compute bistochastic matrix
        bistochastic_matrix = D_r.dot(A).dot(D_c).tocoo()

        return bistochastic_matrix

if __name__ == "__main__":
    normalizer = Normalization()

    bam_file = "./output/alignment/MAP_SORTED.bam"
    fasta_file = "./output/assembly/final_assembly.fasta"
    enzymes = ["EcoRI", "HindIII"]
    output_path = "./output/normalization"

    # Test raw normalization
    raw_output_path = os.path.join(output_path, 'raw')
    raw_cm = normalizer.raw(
        bam_file=bam_file,
        fasta_file=fasta_file,
        enzymes=enzymes,
        output_path=raw_output_path
    )

    if raw_cm is not None:
        print("Raw normalization completed successfully.")
    else:
        print("Raw normalization failed.")

    # Test normcc_contact normalization
    if raw_cm is not None:
        contig_file = os.path.join(raw_output_path, 'contig_info.csv')
        contact_matrix_file = os.path.join(raw_output_path, 'contact_matrix.npz')
        normcc_output_path = os.path.join(output_path, 'normcc')

        norm_result = normalizer.normcc_contact(
            contig_file=contig_file,
            contact_matrix_file=contact_matrix_file,
            output_path=normcc_output_path
        )

        if norm_result is not None:
            print("normCC normalization completed successfully.")
        else:
            print("normCC normalization failed.")
    else:
        print("Skipping normCC normalization due to previous failure.")

    # Test HiCzin normalization
    if raw_cm is not None:
        hiczin_output_path = os.path.join(output_path, 'hiczin')

        hiczin_map = normalizer.hiczin(
            contig_file=contig_file,
            contact_matrix_file=contact_matrix_file,
            output_path=hiczin_output_path,
            min_signal=2  # Adjust as needed
        )

        if hiczin_map is not None:
            print("HiCzin normalization completed successfully.")
        else:
            print("HiCzin normalization failed.")
    else:
        print("Skipping HiCzin normalization due to previous failure.")

    # Test bin3c normalization
    if raw_cm is not None:
        bin3c_output_path = os.path.join(output_path, 'bin3c')

        bin3c_result = normalizer.bin3c(
            contig_file=contig_file,
            contact_matrix_file=contact_matrix_file,
            output_path=bin3c_output_path
        )

        if bin3c_result is not None:
            print("bin3C normalization completed successfully.")
        else:
            print("bin3C normalization failed.")
    else:
        print("Skipping bin3C normalization due to previous failure.")