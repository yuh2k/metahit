import numpy as np
import pandas as pd
import statsmodels.api as sm
import os

def normcc_local(contig_file, output_file=None):
    """
    Perform normcc normalization on the given contig information file.

    Parameters:
    - contig_file: Path to the contig information CSV file.
    - output_file: Path to save the normalization coefficients (optional).

    Returns:
    - norm_result: List of normalization coefficients.
    """
    print(f"Using normcc normalization on contig file: {contig_file}")
    try:
        # Read the contig information
        names = ['contig_name', 'site', 'length', 'covcc', 'signal']
        df = pd.read_csv(
            contig_file,
            header=None,
            names=names,
            dtype={'site': float, 'length': float, 'covcc': float, 'signal': float}
        )

        # 检查数据量
        print(f"Data shape: {df.shape}")
        print(df.head())

        # 检查缺失值
        if df[['site', 'length', 'covcc', 'signal']].isnull().values.any():
            print("Warning: Missing values found.")
            df = df.dropna(subset=['site', 'length', 'covcc', 'signal'])

        # 检查数据量是否足够
        if df.shape[0] <= 4:
            print("Not enough data points to fit the model.")
            return None

        # Apply logarithmic transformation
        epsilon = 1e-6
        df['sample_site'] = np.log(df['site'] + epsilon)
        df['sample_len'] = np.log(df['length'] + epsilon)
        df['sample_covcc'] = np.log(df['covcc'] + epsilon)

        # Define the explanatory and response variables
        exog = df[['sample_site', 'sample_len', 'sample_covcc']]
        endog = df['signal']

        # Add a constant term for the intercept
        exog = sm.add_constant(exog)

        # Fit the Negative Binomial GLM model
        glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial())
        res = glm_nb.fit()

        # Estimate alpha
        alpha = res.scale
        print(f"Estimated alpha: {alpha}")

        # Refit the model with estimated alpha
        glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=alpha))
        res = glm_nb.fit()

        # Get the normalization coefficients
        norm_result = res.params.to_list()

        # Save the normalization coefficients to a file if output_file is specified
        if output_file:
            # 确保输出目录存在
            output_dir = os.path.dirname(output_file)
            os.makedirs(output_dir, exist_ok=True)

            with open(output_file, 'w') as f:
                f.write(','.join(map(str, norm_result)))
                f.write('\n')
            print(f"normcc normalization coefficients saved to {output_file}")

        return norm_result
    except Exception as e:
        print(f"Error during normcc normalization: {e}")
        return None