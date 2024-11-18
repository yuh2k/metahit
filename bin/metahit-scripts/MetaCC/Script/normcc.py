import numpy as np
import pandas as pd
import statsmodels.api as sm

def normcc(contig_file):
    # Define the column names
    names = ['contig_name', 'site', 'length', 'covcc', 'signal']
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(contig_file, header=None, names=names)
    
    # Ensure columns are numeric, coercing any invalid values to NaN
    df['site'] = pd.to_numeric(df['site'], errors='coerce')
    df['length'] = pd.to_numeric(df['length'], errors='coerce')
    df['covcc'] = pd.to_numeric(df['covcc'], errors='coerce')
    df['signal'] = pd.to_numeric(df['signal'], errors='coerce')
    
    # Drop rows with NaN values after conversion
    df.dropna(subset=['site', 'length', 'covcc', 'signal'], inplace=True)
    
    # Add a small constant to avoid log(0) issues
    df['sample_site'] = np.log(df['site'] + 1e-9)
    df['sample_len'] = np.log(df['length'] + 1e-9)
    df['sample_covcc'] = np.log(df['covcc'] + 1e-9)
    
    # Replace any infinite values with NaN
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # Drop rows with NaN or zero values in the signal column
    df.dropna(subset=['sample_site', 'sample_len', 'sample_covcc', 'signal'], inplace=True)
    df = df[df['signal'] > 0]
    
    # Check if DataFrame is empty after preprocessing
    if df.empty:
        raise ValueError("Dataframe is empty after preprocessing or no valid rows with non-zero signal found")

    # Define exogenous (independent) and endogenous (dependent) variables
    exog = df[['sample_site', 'sample_len', 'sample_covcc']]
    endog = df['signal']
    
    # Add a constant term for the regression model
    exog = sm.add_constant(exog)
    
    # Fit a GLM model with Negative Binomial distribution
    try:
        glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
        res = glm_nb.fit(method="lbfgs")
        norm_result = res.params.tolist()
        return norm_result
    except Exception as e:
        print(f"An error occurred during model fitting: {e}")
        return None
