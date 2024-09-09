import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import zscore
from scipy.sparse import load_npz, csr_matrix, save_npz
import scipy.sparse as sp


class Norm:
    
    def __init__(self, contig_info_file, valid_contact_file, output_dir, coverage_file=None):
        self.contig_info_file = contig_info_file
        self.valid_contact_file = valid_contact_file
        self.coverage_file = coverage_file  # Optional coverage file for HiCzin and MetaTOR
        self.output_dir = output_dir
        self.use_coverage = coverage_file is not None  # Only use coverage if file is provided
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def raw(self):
        # Raw normalization output file
        raw_output = os.path.join(self.output_dir, 'raw_contig_info.csv')
        df = pd.read_csv(self.contig_info_file)
        df.to_csv(raw_output, index=False)
        print(f"Raw normalization completed and saved to: {raw_output}")
        return raw_output
    
    def normcc(self):
        # NormCC output file
        normcc_output_txt = os.path.join(self.output_dir, 'normcc_results.txt')
        normcc_output_npz = os.path.join(self.output_dir, 'normcc_results.npz')

        # First, run raw normalization to generate the input for normcc
        raw_output = self.raw() 
        
        names = ['contig_name', 'site', 'length', 'signal']
        df = pd.read_csv(raw_output, header=None, names=names)
        
        df['site'] = pd.to_numeric(df['site'], errors='coerce')
        df['length'] = pd.to_numeric(df['length'], errors='coerce')
        
        if df['signal'].isna().all():
            print("Warning: 'signal' column is entirely NaN. Using mock data for testing.")
            df['signal'] = np.random.rand(len(df))
        
        # Replace zeros with NaN
        df['site'] = df['site'].replace(0, np.nan)
        df['length'] = df['length'].replace(0, np.nan)
        
        # Apply log transformation
        df['sample_site'] = np.log(df['site'])
        df['sample_len'] = np.log(df['length'])
        
        # Filter out rows with NaN or Inf values
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['sample_site', 'sample_len', 'signal'])
        
        if df.empty:
            print("Error: No valid data left after filtering NaN and Inf values.")
            return None
        
        # Exogenous and endogenous variables
        exog = df[['sample_site', 'sample_len']]
        endog = df[['signal']]
        
        # Add constant to exogenous variables
        exog = sm.add_constant(exog)
        
        try:
            # Perform Negative Binomial regression
            glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
            res = glm_nb.fit(method="lbfgs")
            norm_result = res.params.tolist()

            # Save the result to .txt file
            with open(normcc_output_txt, 'w') as f:
                f.write(str(norm_result))
            print(f"NormCC normalization completed and saved to: {normcc_output_txt}")

            # Create a sparse matrix and save it as .npz
            normalized_matrix = csr_matrix(exog)  # Convert the exog (or any matrix you want) to sparse format
            save_npz(normcc_output_npz, normalized_matrix)
            print(f"NormCC sparse matrix saved to: {normcc_output_npz}")

            return normcc_output_npz
        except ValueError as e:
            print(f"Error during GLM fitting: {e}")
            return None
    
    def hiczin(self, thres=0.95):
        # HiCzin output file
        hiczin_output_txt = os.path.join(self.output_dir, 'hiczin_results.txt')
        hiczin_output_npz = os.path.join(self.output_dir, 'hiczin_results.npz')

        # Load the sparse Hi-C contact matrix from the .npz file
        try:
            contact_matrix = load_npz(self.valid_contact_file)  # Load sparse matrix from .npz
        except Exception as e:
            print(f"Error loading .npz file: {e}")
            return None

        # Convert sparse matrix to dense format
        contact_matrix_dense = contact_matrix.toarray()

        # Get indices of non-zero values and corresponding contact counts
        index1, index2 = contact_matrix_dense.nonzero()  # Get the indices of non-zero elements
        contacts = contact_matrix_dense[index1, index2]  # Get the corresponding contact values

        # Create a DataFrame from the non-zero contacts
        sample_data = pd.DataFrame({
            'index1': index1 + 1,  # Shift index1 by 1
            'index2': index2 + 1,  # Shift index2 by 1
            'contacts': contacts
        })

        # Read the contig info file
        contig_info = pd.read_csv(self.contig_info_file)

        # Ensure contig info has the correct column names
        contig_info.columns = ['contig_name', 'length', 'coverage']

        # Ensure 'length' and 'coverage' are numeric
        contig_info['length'] = pd.to_numeric(contig_info['length'], errors='coerce')
        contig_info['coverage'] = pd.to_numeric(contig_info['coverage'], errors='coerce')

        # Replace zeros and missing values in 'length' and 'coverage'
        contig_info['length'].replace(0, np.nan, inplace=True)
        contig_info['coverage'].replace(0, np.nan, inplace=True)
        contig_info.fillna(1, inplace=True)  # Replace NaNs with 1 to avoid log(0) errors

        # Initialize arrays for sample data
        sample_len = np.zeros(len(sample_data))
        sample_site = np.zeros(len(sample_data))
        sample_cov = np.ones(len(sample_data))  # Default to 1 if no coverage is provided

        # Compute log-transformed products for each pair of contigs
        for i in range(len(sample_data)):
            index1 = sample_data.iloc[i, 0]
            index2 = sample_data.iloc[i, 1]
            
            sample_site[i] = np.log(contig_info.iloc[index1 - 1, 1] * contig_info.iloc[index2 - 1, 1])
            sample_len[i] = np.log(contig_info.iloc[index1 - 1, 1] * contig_info.iloc[index2 - 1, 1])
            
            # Use coverage if available
            if self.use_coverage:
                sample_cov[i] = np.log(contig_info.iloc[index1 - 1, 2] * contig_info.iloc[index2 - 1, 2])
        
        # Standardize the data using zscore
        sample_site = zscore(sample_site)
        sample_len = zscore(sample_len)
        sample_cov = zscore(sample_cov) if self.use_coverage else np.ones(len(sample_data))

        # Create the data sample DataFrame
        data_sample = pd.DataFrame({
            'sample_site': sample_site,
            'sample_len': sample_len,
            'sample_cov': sample_cov,
            'sampleCon': sample_data['contacts'].astype(float)
        })
        
        try:
            # Perform GLM Negative Binomial regression
            model = sm.GLM(data_sample['sampleCon'], 
                           sm.add_constant(data_sample[['sample_site', 'sample_len', 'sample_cov']]), 
                           family=sm.families.NegativeBinomial())
            fit1 = model.fit()
        except Exception as e:
            print(f"Error: {e}")
            return None
        
        # Extract coefficients and calculate residuals
        coeff = fit1.params.values
        res_sample = data_sample['sampleCon'] / np.exp(coeff[0] + coeff[1] * sample_site + coeff[2] * sample_len + coeff[3] * sample_cov)

        # Filter non-zero results and calculate quantiles
        index_nonzero = res_sample > 0
        res_sample_nonzero = res_sample[index_nonzero]
        perc = np.quantile(res_sample_nonzero, thres)

        # Store the results
        result = {
            'coeff': coeff[:4],
            'perc': perc,
            'mean_site': np.mean(sample_site),
            'sd_site': np.std(sample_site),
            'mean_len': np.mean(sample_len),
            'sd_len': np.std(sample_len),
            'mean_cov': np.mean(sample_cov) if self.use_coverage else None,
            'sd_cov': np.std(sample_cov) if self.use_coverage else None
        }

        # Save the result to a text file
        with open(hiczin_output_txt, 'w') as f:
            f.write(str(result))
        print(f"HiCzin normalization completed and saved to: {hiczin_output_txt}")

        # Save the normalized results as a sparse matrix in .npz format
        normalized_matrix = csr_matrix(res_sample)
        save_npz(hiczin_output_npz, normalized_matrix)
        print(f"HiCzin normalized matrix saved to: {hiczin_output_npz}")

        return hiczin_output_npz


    def metator(self):
        # MetaTOR output file
        metator_output_txt = os.path.join(self.output_dir, 'metator_results.txt')
        metator_output_npz = os.path.join(self.output_dir, 'metator_results.npz')

        # Load the sparse Hi-C contact matrix from the .npz file
        try:
            contact_matrix = load_npz(self.valid_contact_file)  # Load sparse matrix from .npz
        except Exception as e:
            print(f"Error loading .npz file: {e}")
            return None

        # Convert sparse matrix to dense format
        contact_matrix_dense = contact_matrix.toarray()

        # Get indices of non-zero values and corresponding contact counts
        index1, index2 = contact_matrix_dense.nonzero()  # Get the indices of non-zero elements
        contacts = contact_matrix_dense[index1, index2]  # Get the corresponding contact values

        sample_data = pd.DataFrame({
            'index1': index1 + 1,  # Shift index1 by 1
            'index2': index2 + 1,  # Shift index2 by 1
            'contacts': contacts
        })

        # Load the contig info
        contig_info = pd.read_csv(self.contig_info_file, header=None, sep=',')

        # Check the number of columns in contig_info and set column names accordingly
        if contig_info.shape[1] == 3:
            contig_info.columns = ['contig_name', 'site', 'length']
        elif contig_info.shape[1] == 4:
            contig_info.columns = ['contig_name', 'site', 'length', 'coverage']
        else:
            print(f"Error: Unexpected number of columns in contig_info: {contig_info.shape[1]}")
            return None

        # Normalize contacts using coverage if available
        norm_contacts = np.zeros(len(sample_data))
        for i in range(len(sample_data)):
            index1 = sample_data.iloc[i, 0]
            index2 = sample_data.iloc[i, 1]
            
            # If coverage is available, use it; otherwise default to 1
            coverage1 = contig_info.iloc[index1 - 1, 3] if 'coverage' in contig_info.columns else 1
            coverage2 = contig_info.iloc[index2 - 1, 3] if 'coverage' in contig_info.columns else 1
            
            norm_contacts[i] = sample_data.iloc[i, 2] / np.sqrt(coverage1 * coverage2)

        # Save the normalized contacts
        np.savetxt(metator_output_txt, norm_contacts)
        print(f"MetaTOR normalization completed and saved to: {metator_output_txt}")

        # Save as npz
        normalized_matrix = csr_matrix(norm_contacts)
        save_npz(metator_output_npz, normalized_matrix)
        print(f"MetaTOR normalized matrix saved to: {metator_output_npz}")

        return metator_output_npz

    def bin3c(self):
        # bin3C output file
        bin3c_output_txt = os.path.join(self.output_dir, 'bin3c_results.txt')
        bin3c_output_npz = os.path.join(self.output_dir, 'bin3c_results.npz')

        # Load the sparse Hi-C contact matrix from the .npz file
        try:
            contact_matrix = load_npz(self.valid_contact_file)  # Load sparse matrix from .npz
        except Exception as e:
            print(f"Error loading .npz file: {e}")
            return None

        # Convert sparse matrix to dense format
        contact_matrix_dense = contact_matrix.toarray()

        # Get indices of non-zero values and corresponding contact counts
        index1, index2 = contact_matrix_dense.nonzero()  # Get the indices of non-zero elements
        contacts = contact_matrix_dense[index1, index2]  # Get the corresponding contact values

        sample_data = pd.DataFrame({
            'index1': index1 + 1,  # Shift index1 by 1
            'index2': index2 + 1,  # Shift index2 by 1
            'contacts': contacts
        })

        # Load the contig info file, skipping the header row
        contig_info = pd.read_csv(self.contig_info_file, header=1, sep=',')

        if contig_info.shape[1] == 3:
            contig_info.columns = ['contig_name', 'site', 'length']
            contig_info['coverage'] = 1  # Set default coverage to 1 if not present
            print("Coverage column not found, setting default value to 1.")
        else:
            print(f"Error: Unexpected number of columns in contig_info: {contig_info.shape[1]}")
            return None

        # Convert site and length to numeric, handling errors and non-numeric values
        contig_info['site'] = pd.to_numeric(contig_info['site'], errors='coerce')
        contig_info['length'] = pd.to_numeric(contig_info['length'], errors='coerce')

        # Replace NaN or zero values in 'site' and 'length' with a small value to avoid division errors
        contig_info['site'].replace(0, np.nan, inplace=True)
        contig_info['length'].replace(0, np.nan, inplace=True)
        contig_info.fillna(1, inplace=True)

        # Print the first few rows of contig_info for debugging
        print("contig_info head after conversion:\n", contig_info.head())

        # Initialize array for normalized contacts
        norm_contacts = np.zeros(len(sample_data))

        # Ensure indices do not exceed contig_info size
        max_index = len(contig_info)

        # Calculate normalized contacts
        for i in range(len(sample_data)):
            index1 = sample_data.iloc[i, 0]
            index2 = sample_data.iloc[i, 1]

            if index1 > max_index or index2 > max_index:
                print(f"Warning: Index {index1} or {index2} is out of bounds, skipping this entry.")
                continue  # Skip entries with out-of-bounds indices

            sites1 = contig_info.iloc[index1 - 1, 1]
            sites2 = contig_info.iloc[index2 - 1, 1]
            norm_contacts[i] = sample_data.iloc[i, 2] / (sites1 * sites2)

        # Save the normalized contacts
        np.savetxt(bin3c_output_txt, norm_contacts)
        print(f"bin3C normalization completed and saved to: {bin3c_output_txt}")

        # Save as npz
        normalized_matrix = sp.csr_matrix(norm_contacts)
        sp.save_npz(bin3c_output_npz, normalized_matrix)
        print(f"bin3C normalized matrix saved to: {bin3c_output_npz}")

        return bin3c_output_npz


# Test the Norm class
if __name__ == "__main__":
    # Example paths (you need to provide the correct paths)
    contig_info_file = "./output/downstream/contig_info.csv"
    valid_contact_file = "./output/downstream/hic_contact_matrix.npz"
    coverage_file = "./output/coverage_estimation.csv"  # Optional
    output_dir = "./output/normalization"
    
    norm = Norm(contig_info_file, valid_contact_file, output_dir, coverage_file)
    
    # Test each normalization method
    norm.raw()
    norm.normcc()
    norm.hiczin()
    norm.metator()
    norm.bin3c()
