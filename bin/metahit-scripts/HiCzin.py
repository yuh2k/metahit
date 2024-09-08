import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import zscore

def HiCzin(contig_info_file, valid_contact_file, thres):
    # Read the valid contact file and contig info file
    sample_data = pd.read_csv(valid_contact_file, header=None, sep=',')
    sample_data.columns = ['index1', 'index2', 'contacts']
    contig_info = pd.read_csv(contig_info_file, header=None, sep=',')
    
    # Adjust the indexes by adding 1
    sample_data['index1'] = sample_data['index1'] + 1
    sample_data['index2'] = sample_data['index2'] + 1
    
    if contig_info.shape[1] == 4:
        # Set the column names for contig info
        contig_info.columns = ['contig_name', 'site', 'length', 'coverage']
        
        # Replace 0 values for site with 1 and replace 0 values for coverage with the minimum non-zero coverage
        contig_info.loc[contig_info['site'] == 0, 'site'] = 1
        contig_info.loc[contig_info['coverage'] == 0, 'coverage'] = contig_info[contig_info['coverage'] != 0]['coverage'].min()
        
        # Initialize arrays for site, length, and coverage
        sample_len = np.zeros(len(sample_data))
        sample_site = np.zeros(len(sample_data))
        sample_cov = np.zeros(len(sample_data))
        
        # Compute the log-transformed product of site, length, and coverage for each pair of contigs
        for i in range(len(sample_data)):
            index1 = sample_data.iloc[i, 0]
            index2 = sample_data.iloc[i, 1]
            
            sample_site[i] = np.log(contig_info.iloc[index1 - 1, 1] * contig_info.iloc[index2 - 1, 1])
            sample_len[i] = np.log(contig_info.iloc[index1 - 1, 2] * contig_info.iloc[index2 - 1, 2])
            sample_cov[i] = np.log(contig_info.iloc[index1 - 1, 3] * contig_info.iloc[index2 - 1, 3])
        
        sampleCon = sample_data['contacts'].astype(float).values
        
        sample_site = zscore(sample_site)
        sample_len = zscore(sample_len)
        sample_cov = zscore(sample_cov)
        
        data_sample = pd.DataFrame({
            'sample_site': sample_site,
            'sample_len': sample_len,
            'sample_cov': sample_cov,
            'sampleCon': sampleCon
        })
        
        try:
            model = sm.GLM(data_sample['sampleCon'], 
                           sm.add_constant(data_sample[['sample_site', 'sample_len', 'sample_cov']]), 
                           family=sm.families.NegativeBinomial())
            fit1 = model.fit()
        except Exception as e:
            print(f"Error: {e}")
            return None
        
        coeff = fit1.params.values

        res_sample = sampleCon / np.exp(coeff[0] + coeff[1] * sample_site + coeff[2] * sample_len + coeff[3] * sample_cov)

        index_nonzero = res_sample > 0
        res_sample_nonzero = res_sample[index_nonzero]

        perc = np.quantile(res_sample_nonzero, thres)
    
        result = {
            'coeff': coeff[:4],
            'perc': perc,
            'mean_site': np.mean(sample_site),
            'sd_site': np.std(sample_site),
            'mean_len': np.mean(sample_len),
            'sd_len': np.std(sample_len),
            'mean_cov': np.mean(sample_cov),
            'sd_cov': np.std(sample_cov)
        }
        
        return result
