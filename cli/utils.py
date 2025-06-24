"""
    Utility functions for the CLI.
    This module contains functions to handle the command line interface (CLI) for the application.
    It includes functions to parse command line arguments, display help messages, and handle errors.
"""
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


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

def plot_kl_distance(df):
    """
    Plot stacked bar chart of KL distances by number of clusters.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with KL distance data where:
        - 'cluster' column contains the number of clusters
        - 'group1' through 'groupN' columns contain individual KL distances
        - 'total' column contains the total KL distance
    
    Returns:
    --------
    fig, ax : matplotlib figure and axes
        The plotted figure and axes objects
    """
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract the number of clusters and group columns
    cluster_nums = df['cluster'].values
    group_cols = [col for col in df.columns if col.startswith('group')]
    
    # Create bottom positions for stacking
    bottoms = np.zeros(len(df))
    
    # Plot each group as a segment of the stacked bar
    for i, col in enumerate(group_cols):
        # Skip if the column contains only zeros
        if df[col].sum() == 0:
            continue
            
        values = df[col].values
        # Only plot non-zero values
        mask = values > 0
        
        if mask.any():
            bars = ax.bar(
                cluster_nums[mask], 
                values[mask], 
                bottom=bottoms[mask],
                width=0.6, 
                color='black', 
                edgecolor='white'
            )
            
            # Add group numbers as text in the middle of each bar segment
            for j, (x, y, height) in enumerate(zip(cluster_nums[mask], bottoms[mask], values[mask])):
                if height > 0.5:  # Only add text if the segment is large enough
                    ax.text(
                        x, 
                        y + height/2, 
                        f"{i+1}", 
                        ha='center', 
                        va='center', 
                        color='white', 
                        fontsize=12
                    )
            
            # Update bottoms for next stack
            bottoms[mask] += values[mask]
    
    # Set labels and title
    ax.set_xlabel('Number of clusters', fontsize=12)
    ax.set_ylabel('Kullback Leibler distance', fontsize=12)
    
    # Set y-axis limits
    ax.set_ylim(0, max(df['total']) * 1.05)
    
    # Set x-axis ticks to be the cluster numbers
    ax.set_xticks(cluster_nums)
    
    # Return the figure and axis
    return fig, ax

# fig,ax =  plot_kl_distance(read_KLD_file("/Users/sanjay/Monash/Master_thesis/lab_work/Li_Lab/cluster/HLA-PepClust/data/D90_HLA_3844874/images/gibbs.KLDvsClusters.tab"))
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_kld_pcc_distance(df, correlation_df=None):
    """
    Plot stacked bar chart of KL distances by number of clusters,
    with HLA and correlation information inside each bar segment.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with KL distance data where:
        - 'cluster' column contains the number of clusters
        - 'group1' through 'groupN' columns contain individual KL distances
        - 'total' column contains the total KL distance
    
    correlation_df : pandas.DataFrame, optional
        DataFrame with correlation data where:
        - 'Cluster' column contains cluster identifiers (e.g., '1of3')
        - 'HLA' column contains HLA identifiers
        - 'Correlation' column contains correlation values
        - 'KLD' column contains KL distance values
    
    Returns:
    --------
    fig, ax : matplotlib figure and axes
        The plotted figure and axes objects
    """
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Extract the number of clusters and group columns
    cluster_nums = df['cluster'].values
    group_cols = [col for col in df.columns if col.startswith('group')]
    
    # Create a dictionary to store best correlation data for each group in each cluster
    best_corr_by_group = {}
    
    if correlation_df is not None:
        # Prepare correlation lookup
        for cluster_num in cluster_nums:
            best_corr_by_group[cluster_num] = {}
            
            for group_num in range(1, cluster_num + 1):
                # Format for the cluster identifier - e.g., "1of3" for group 1 in 3-cluster solution
                cluster_pattern = f"{group_num}of{cluster_num}"
                
                # Find this group in the correlation data
                group_data = correlation_df[correlation_df['Cluster'] == cluster_pattern]
                print(group_data)
                if not group_data.empty:
                    # Get the row with the highest correlation
                    best_row = group_data.loc[group_data['Correlation'].idxmax()]
                    best_corr_by_group[cluster_num][group_num] = {
                        'hla': best_row['HLA'],
                        'correlation': best_row['Correlation'],
                        'kld': best_row['KLD']
                    }
    
    # Create bottom positions for stacking
    bottoms = np.zeros(len(df))
    
    # Plot each group as a segment of the stacked bar
    for i, col in enumerate(group_cols):
        # Skip if the column contains only zeros
        if df[col].sum() == 0:
            continue
            
        values = df[col].values
        # Only plot non-zero values
        mask = values > 0
        
        if mask.any():
            bars = ax.bar(
                cluster_nums[mask], 
                values[mask], 
                bottom=bottoms[mask],
                width=0.6, 
                color='black', 
                edgecolor='white'
            )
            
            # Add group numbers and correlation info in the middle of each bar segment
            for j, (x, y, height) in enumerate(zip(cluster_nums[mask], bottoms[mask], values[mask])):
                # Get the cluster number and group number
                cluster_num = x
                group_num = i + 1
                
                # Only add text if the segment is large enough
                if height > 0.5:
                    # First add the group number
                    ax.text(
                        x, 
                        y + height/2,
                        f"{group_num}", 
                        ha='center', 
                        va='center', 
                        color='white', 
                        fontsize=14,
                        fontweight='bold'
                    )
                    
                    # Add HLA and correlation if available
                    if correlation_df is not None and cluster_num in best_corr_by_group and group_num in best_corr_by_group[cluster_num]:
                        # Get the correlation data for this group
                        corr_data = best_corr_by_group[cluster_num][group_num]
                        # print(corr_data)
                        # Add text at the top of the current bar segment
                        if height > 1.2:  # Only add if there's enough room
                            ax.text(
                                x,
                                y + height - 0.2,
                                f"{str(corr_data['hla']).replace("_","*")}\npcc = {corr_data['correlation']:.2f}", #kld = {corr_data['kld']:.2f}",
                                ha='center',
                                va='top',
                                color='white',
                                fontsize=8,
                                fontweight='bold'
                            )
            
            # Update bottoms for next stack
            bottoms[mask] += values[mask]
    
    # Add the highest correlation HLA at the top of each complete bar
    # if correlation_df is not None:
    #     for cluster_num in cluster_nums:
    #         # Format for the cluster identifier - e.g., "4of4" for the complete 4-cluster solution
    #         cluster_pattern = f"{cluster_num}of{cluster_num}"
            
    #         # Find this cluster in the correlation data
    #         cluster_data = correlation_df[correlation_df['Cluster'] == cluster_pattern]
            
    #         if not cluster_data.empty:
    #             # Get the row with the highest correlation
    #             best_row = cluster_data.loc[cluster_data['Correlation'].idxmax()]
                
    #             # Get the height of the corresponding bar
    #             bar_height = df.loc[df['cluster'] == cluster_num, 'total'].values[0]
                
    #             # Add text above the bar
    #             ax.text(
    #                 cluster_num,
    #                 bar_height * 1.02,
    #                 f"{best_row['HLA']}\npcc = {best_row['Correlation']:.2f}",
    #                 ha='center',
    #                 va='bottom',
    #                 fontsize=11,
    #                 fontweight='bold'
    #             )
    
    # Set labels and title
    ax.set_xlabel('Number of clusters', fontsize=12)
    ax.set_ylabel('Kullback Leibler distance', fontsize=12)
    
    # Set y-axis limits with extra space for correlation text
    max_height = max(df['total'])
    ax.set_ylim(0, max_height * 1.15)
    
    # Set x-axis ticks to be the cluster numbers
    ax.set_xticks(cluster_nums)
    
    # Add a tight layout
    plt.tight_layout()
    
    # Return the figure and axis
    return fig, ax

