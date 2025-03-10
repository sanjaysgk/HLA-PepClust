"""
Utility functions for HLA-PEPCLUST.

This module contains utility functions used across the application.
"""

import os
import re
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any


def generate_unique_random_ids(count: int) -> List[int]:
    """
    Generate a list of unique 6-digit random IDs.

    Args:
        count (int): Number of unique random IDs to generate.
        
    Returns:
        List[int]: List of unique 6-digit random IDs.
        
    Raises:
        ValueError: If count is larger than the range of available IDs.
    """
    start = 100000
    end = 999999
    
    if count > (end - start + 1):
        raise ValueError("Count is larger than the range of unique IDs available.")

    return np.random.choice(range(start, end + 1), size=count, replace=False).tolist()


def ensure_directory_exists(path: str) -> str:
    """
    Create directory if it doesn't exist.
    
    Args:
        path (str): Path to the directory to create.
        
    Returns:
        str: Path to the created directory.
    """
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def format_hla_name(hla_name: str, for_db: bool = False) -> str:
    """
    Format HLA name for consistent use.
    
    Args:
        hla_name (str): The HLA name to format.
        for_db (bool): If True, formats for database comparison.
        
    Returns:
        str: Formatted HLA name.
    """
    formatted = hla_name.replace("HLA_", "").replace("*", "").replace("txt", "").replace(".", "")
    return formatted


def parse_hla_list(hla_list: str) -> List[str]:
    """
    Parse a comma-separated string of HLA types into a list.
    
    Args:
        hla_list (str): Comma-separated string of HLA types.
        
    Returns:
        List[str]: List of formatted HLA types.
    """
    if not hla_list:
        return []
        
    hla_types = hla_list.split(",")
    return [format_hla_name(hla) for hla in hla_types]


def format_input_gibbs(gibbs_matrix: str) -> pd.DataFrame:
    """
    Format the Gibbs output matrix.

    Args:
        gibbs_matrix (str): Path to the Gibbs matrix file.
        
    Returns:
        pd.DataFrame: Processed DataFrame.
    """
    df = pd.read_csv(gibbs_matrix)
    amino_acids = df.iloc[0, 0].split()
    df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"^\d+\s\w\s", "", regex=True)
    
    new_df = pd.DataFrame(
        df.iloc[1:, 0].str.split(expand=True).values, 
        columns=amino_acids
    )
    
    new_df.reset_index(drop=True, inplace=True)
    new_df = new_df.apply(pd.to_numeric, errors='coerce')
    
    return new_df


def amino_acid_order_identical(df1: pd.DataFrame, df2: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Ensure amino acid column order is identical between two dataframes.
    
    Args:
        df1 (pd.DataFrame): First DataFrame with amino acid columns.
        df2 (pd.DataFrame): Second DataFrame with amino acid columns.
        
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Tuple of DataFrames with aligned columns.
    """
    if list(df1.columns) != list(df2.columns):
        # Reorder the second DataFrame's columns to match the first
        df2 = df2[df1.columns]
    return df1, df2


def find_file_by_pattern(directory: str, pattern: str, verbose: bool = False) -> str:
    """
    Find a file in a directory that matches the given pattern.
    
    Args:
        directory (str): Directory to search in.
        pattern (str): Pattern to match in filenames.
        verbose (bool): If True, print debugging information.
        
    Returns:
        str: Full path to the matched file, or None if not found.
    """
    if not os.path.exists(directory):
        if verbose:
            print(f"Directory does not exist: {directory}")
        return None
    
    if verbose:
        print(f"Searching for '{pattern}' in {directory}")
        print(f"Directory contents: {os.listdir(directory)}")
    
    for filename in os.listdir(directory):
        if pattern in filename:
            return os.path.join(directory, filename)
    
    if verbose:
        print(f"No file matching '{pattern}' found in {directory}")
    return None


def create_output_directories(base_path: str) -> Dict[str, str]:
    """
    Create standard output directory structure for results.
    
    Args:
        base_path (str): Base path for output directories.
        
    Returns:
        Dict[str, str]: Dictionary of created directory paths.
    """
    dirs = {
        "root": ensure_directory_exists(base_path),
        "cluster_img": ensure_directory_exists(os.path.join(base_path, "cluster-img")),
        "allotypes_img": ensure_directory_exists(os.path.join(base_path, "allotypes-img")),
        "corr_data": ensure_directory_exists(os.path.join(base_path, "corr-data"))
    }
    
    return dirs