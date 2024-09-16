import numpy as np
import pandas as pd
import statsmodels.api as sm
from norm.normcc import *
from norm.bin3c import ContactMap  
from norm.hiczin import *          
import scipy.sparse as scisp
from scipy.sparse import csr_matrix, save_npz

class Normalization:
    def __init__(self):
        pass

    # Method for raw normalization (pass-through method)
    def raw(self, file, output_file):
        print(f"Using raw normalization on file: {file}")
        try:
            df = pd.read_csv(file)
            df.to_csv(output_file, index=False)
            print(f"Raw normalization result saved to {output_file}")
            return df
        except Exception as e:
            print(f"Error during raw normalization: {e}")
            return None

    # Method for normcc normalization
    def normcc(self, contig_file, output_file):
        print(f"Applying NormCC normalization on {contig_file}")
        try:
            names = ['contig_name', 'site', 'length', 'covcc', 'signal']
            df = pd.read_csv(contig_file, header=None, names=names)

            # Ensure 'site', 'length', and 'covcc' are numeric and handle zero/negative values
            df['site'] = pd.to_numeric(df['site'], errors='coerce').replace(0, 1e-6)
            df['length'] = pd.to_numeric(df['length'], errors='coerce').replace(0, 1e-6)
            df['covcc'] = pd.to_numeric(df['covcc'], errors='coerce').replace(0, 1e-6)

            # Check for empty arrays and log them
            if df['site'].isnull().all() or df['length'].isnull().all() or df['covcc'].isnull().all():
                raise ValueError("One or more columns in the dataframe are empty after cleaning!")

            # Apply log transformation
            df['sample_site'] = np.log(df['site'])
            df['sample_len'] = np.log(df['length'])
            df['sample_covcc'] = np.log(df['covcc'])

            # Drop rows with NaN values created during the log transformation
            df.dropna(subset=['sample_site', 'sample_len', 'sample_covcc'], inplace=True)

            # Prepare the model
            exog = df[['sample_site', 'sample_len', 'sample_covcc']]
            endog = df[['signal']]
            exog = sm.add_constant(exog)

            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit(method="lbfgs")
            norm_result = res.params.tolist()

            # Save the result
            df['norm_signal'] = res.predict(exog)
            df.to_csv(output_file, index=False)
            print(f"NormCC normalization result saved to {output_file}")
            return norm_result
        except Exception as e:
            print(f"Error during NormCC normalization: {e}")
            return None

    # Method for bin3c normalization
    def bin3c(self, bam_file, seq_file, enzymes, bin_size=1000, coverage_file=None, output_file="bin3c_output.npz"):
        print(f"Applying bin3c normalization on BAM: {bam_file}, using sequence file: {seq_file}")
        try:
            # Ensure enzymes is a list
            if enzymes is None or not isinstance(enzymes, list):
                raise ValueError("enzyme_names must be a collection of names")
            
            contact_map = ContactMap(bam_file, enzymes, seq_file, min_insert=500, bin_size=bin_size)
            contact_map.prepare_seq_map()

            processed_map = contact_map.processed_map  # Assuming this is a sparse matrix or can be converted to one

            # Check if processed_map is a sparse matrix, if not convert it
            if not scisp.issparse(processed_map):
                processed_map = csr_matrix(processed_map)
            
            # Save the sparse matrix
            save_npz(output_file, processed_map)
            
            print(f"Bin3C normalization result saved to {output_file}")
            return processed_map
        except Exception as e:
            print(f"Error during bin3c normalization: {e}")
            return None

    # Method for hiczin normalization
    def hiczin(self, bam_file, seq_file, enzymes, bin_size=1000, tip_size=500, output_file="hiczin_output.npz"):
        print(f"Applying hicZin normalization on BAM: {bam_file}, using sequence file: {seq_file}")
        try:
            # Ensure enzymes is a list
            if enzymes is None or not isinstance(enzymes, list):
                raise ValueError("enzyme_names must be a collection of names")

            contact_map = ContactMap(bam_file, enzymes, seq_file, min_insert=500, bin_size=bin_size, tip_size=tip_size)
            contact_map.prepare_seq_map()

            processed_map = contact_map.processed_map  # Assuming this is a sparse matrix or can be converted to one

            # Check if processed_map is a sparse matrix, if not convert it
            if not scisp.issparse(processed_map):
                processed_map = csr_matrix(processed_map)
            
            # Save the sparse matrix
            save_npz(output_file, processed_map)
            
            print(f"HicZin normalization result saved to {output_file}")
            return processed_map
        except Exception as e:
            print(f"Error during HicZin normalization: {e}")
            return None


# Example of how to use the Normalization class
if __name__ == "__main__":
    normalizer = Normalization()

    # Example of raw normalization
    # raw_result = normalizer.raw("./output/downstream/contig_info.csv", "./output/downstream/raw_output.csv")
    
    # if raw_result is not None:
    #     normcc_result = normalizer.normcc("./output/downstream/contig_info.csv", "./output/downstream/normcc_output.csv")

    # # Test bin3c and hiczin normalization
    enzymes = ["EcoRI", "HindIII"]  # Example enzyme names
    # bin3c_result = normalizer.bin3c("./output/alignment/MAP_SORTED.bam", "./output/assembly/final_assembly.fasta", enzymes, output_file="./output/bin3c_output.npz")
    hiczin_result = normalizer.hiczin("./output/alignment/MAP_SORTED.bam", "./output/assembly/final_assembly.fasta", enzymes, output_file="./output/hiczin_output.npz")
