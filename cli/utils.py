"""
    Utility functions for the CLI.
    This module contains functions to handle the command line interface (CLI) for the application.
    It includes functions to parse command line arguments, display help messages, and handle errors.
"""
import os
import sys
import pandas as pd


def check_python_version():
    """
    Check if the Python version is 3.8 or higher.
    If not, print an error message and exit the program.
    """
    if sys.version_info < (3, 8):
        print("Error: Python 3.8 or higher is required.")
        sys.exit(1)
    
    if sys.version_info >= (3, 10):
        print("Python version is 3.10 or higher.")
        
def check_environment_variable(var_name):
    """
    Check if a specific environment variable is set.
    If not, print an error message and exit the program.
    
    Args:
        var_name (str): The name of the environment variable to check.
    """
    
    if var_name not in os.environ:
        print(f"Error: Environment variable '{var_name}' is not set.")
        sys.exit(1)
    else:
        print(f"Environment variable '{var_name}' is set to: {os.environ[var_name]}")
        
def read_KLD_file(file_path):
    """
    Read & Prase Gibbs gibbs.KLDvsClusters file.
    Args:
        file_path (str): The path to the Gibbs gibbs.KLDvsClusters file.
    Returns:
        list: A list of tuples containing the cluster number and KLD value.
    """
    try:
        # Initialize an empty dictionary to store data
        data_dict = {'cluster': []}
        
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines[1:]:  # Skip the header line
                parts = line.strip().split('\t')
                cluster_number = int(parts[0])
                kld_values = [float(value) for value in parts[1:]]
                
                # Ensure the dictionary has enough group columns
                for i in range(1, len(kld_values) + 1):
                    group_key = f'group{i}'
                    if group_key not in data_dict:
                        data_dict[group_key] = []
                
                # Populate the row data
                data_dict['cluster'].append(cluster_number)
                for i, kld_value in enumerate(kld_values, start=1):
                    data_dict[f'group{i}'].append(kld_value if kld_value > 0 else 0)
                
                # Fill remaining groups with zeros if necessary
                for j in range(len(kld_values) + 1, len(data_dict) - 1):
                    data_dict[f'group{j}'].append(0)
        
        # Convert the dictionary to a DataFrame
        df = pd.DataFrame(data_dict)
        df.loc[:, 'total'] = df.iloc[:, 1:].sum(axis=1)
        df.reset_index(drop=True, inplace=True)
        # print(df)
        return df
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

# read_KLD_file("/Users/sanjay/Monash/Master_thesis/lab_work/Li_Lab/cluster/HLA-PepClust/data/D90_HLA_3844874/images/gibbs.KLDvsClusters.tab")