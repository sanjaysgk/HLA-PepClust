
"""
Configuration management for HLA-PEPCLUST.

This module handles configuration loading, validation, and default settings.
"""

import os
import json
from typing import Dict, Any
from cli.logger import CONSOLE
from importlib.metadata import distributions

for dist in distributions():
    print(dist.metadata["Name"], dist.metadata["Version"])


# Default configuration if no config file is provided
DEFAULT_CONFIG = {
    "human": {
        "path": "Gibbs_motifs_human",
        "matrix": "matrices",
        "motif": "motif",
        "allotypes": "allotypes/hla_data.csv",
        "ref_data": "data/ref_data"
    },
    "mouse": {
        "path": "Gibbs_motifs_mouse",
        "matrix": "matrices",
        "motif": "motif",
        "allotypes": "allotypes/mhc_data.csv",
        "ref_data": "data/ref_data"
    }
}


def load_config(config_file: str) -> Dict[str, Any]:
    """
    Load configuration from a JSON file, or return default if file doesn't exist.
    
    Args:
        config_file (str): Path to the configuration file.
        
    Returns:
        Dict[str, Any]: Configuration dictionary.
    """
    if not os.path.exists(config_file):
        CONSOLE.log(f"Config file {config_file} does not exist", style="red")
        CONSOLE.log("Using default configuration", style="blue")
        return DEFAULT_CONFIG
    
    try:
        with open(config_file, "r") as f:
            config = json.load(f)
        return config
    except json.JSONDecodeError:
        CONSOLE.log(f"Error parsing config file {config_file}. Using default configuration.", style="red")
        return DEFAULT_CONFIG


def validate_config(config: Dict[str, Any]) -> bool:
    """
    Validate the configuration structure.
    
    Args:
        config (Dict[str, Any]): Configuration dictionary.
        
    Returns:
        bool: True if valid, False otherwise.
    """
    required_species = ["human", "mouse"]
    required_keys = ["path", "matrix", "motif", "allotypes", "ref_data"]
    
    for species in required_species:
        if species not in config:
            CONSOLE.log(f"Missing required species '{species}' in configuration", style="red")
            return False
        
        for key in required_keys:
            if key not in config[species]:
                CONSOLE.log(f"Missing required key '{key}' for species '{species}' in configuration", style="red")
                return False
    
    return True


def verify_paths(config: Dict[str, Any]) -> bool:
    """
    Verify that all paths in the configuration exist.
    
    Args:
        config (Dict[str, Any]): Configuration dictionary.
        
    Returns:
        bool: True if all paths exist, False otherwise.
    """
    for species, species_config in config.items():
        # Check main directory
        species_path = os.path.join(species_config["ref_data"], species_config["path"])
        if not os.path.exists(species_path):
            CONSOLE.log(f"Path {species_path} does not exist. Check the path.", style="red")
            return False
        
        # Check matrix directory
        matrix_path = os.path.join(species_path, species_config["matrix"])
        if not os.path.exists(matrix_path):
            CONSOLE.log(f"Matrix path {matrix_path} does not exist", style="red")
            return False
        
        # Check motif directory
        motif_path = os.path.join(species_path, species_config["motif"])
        if not os.path.exists(motif_path):
            CONSOLE.log(f"Motif path {motif_path} does not exist", style="red")
            return False
        
        # Check allotypes file
        allotypes_path = os.path.join(species_config["ref_data"], species_config["allotypes"])
        if not os.path.exists(allotypes_path):
            CONSOLE.log(f"Allotype path {allotypes_path} does not exist", style="red")
            return False
    
    return True