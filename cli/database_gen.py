"""
Database generation module for HLA-PEPCLUST.

This module handles the generation of reference database files for HLA/MHC allotypes.
"""

import os
import sys
import time
import pandas as pd
from typing import Dict, Any

from cli.logger import CONSOLE
from cli.config import load_config, validate_config, verify_paths


def process_allotype_data(config: Dict[str, Any]) -> None:
    """
    Process allotype data and create database files.
    
    Args:
        config (Dict[str, Any]): Configuration dictionary.
    """
    for species, species_config in config.items():
        # Load allotype data
        allotype_file = os.path.join(species_config["ref_data"], species_config["allotypes"])
        try:
            allotypes = pd.read_csv(allotype_file)
        except Exception as e:
            CONSOLE.log(f"Error reading allotype file {allotype_file}: {str(e)}", style="red")
            sys.exit(1)
            
        # Process based on species
        if species == "human":
            allotypes.rename(columns={"formatted_HLA": "formatted_allotypes", "HLA": "allotypes"}, inplace=True)
            allotypes['motif'] = "HLA_" + allotypes['formatted_allotypes'] + ".png"
        elif species == "mouse":
            allotypes.rename(columns={"formatted_MHC": "formatted_allotypes", "MHC": "allotypes"}, inplace=True)
            allotypes['motif'] = allotypes['formatted_allotypes'] + ".png"
        else:
            CONSOLE.log(f"Species {species} not supported", style="red")
            sys.exit(1)
            
        # Add species and initialize path columns
        allotypes['species'] = species
        allotypes['motif_path'] = ""
        allotypes['matrices_path'] = ""
        
        # Verify and add paths for each motif
        for motif in allotypes['motif']:
            motif_path = os.path.join(
                species_config['ref_data'], 
                species_config['path'], 
                species_config['motif'], 
                motif
            )
            matrix_path = os.path.join(
                species_config['ref_data'],
                species_config['path'],
                species_config['matrix'],
                motif.replace('.png', '.txt')
            )
            
            CONSOLE.log(f"[yellow]{species}[/yellow] checking {motif_path}", style="blue")
            
            if not os.path.exists(motif_path):
                CONSOLE.log(f"Motif {motif} does not exist at {motif_path}", style="red")
                sys.exit(1)
            elif not os.path.exists(matrix_path):
                CONSOLE.log(f"Matrix {motif.replace('.png', '.txt')} does not exist at {matrix_path}", style="red")
                sys.exit(1)
            else:
                # Update paths in dataframe - store absolute paths to ensure consistent access
                motif_path = os.path.join(species_config['path'], 'motif', motif)
                matrix_path = os.path.join(species_config['path'], 'matrices', motif.replace('.png', '.txt'))
                
                # Store these paths - important to use consistent path format
                allotypes.loc[allotypes['motif'] == motif, 'motif_path'] = motif_path
                allotypes.loc[allotypes['motif'] == motif, 'matrices_path'] = matrix_path
                
                CONSOLE.log(f"[green]✓[/green] Registered {motif} with matrix at {matrix_path}", style="dim")
        
        # Save database
        db_path = os.path.join(species_config['ref_data'], f"{species}.db")
        allotypes.to_csv(db_path, index=False)
        CONSOLE.log(f"Database {species}.db created successfully ({db_path})", style="green")


def generate_database(config_file: str) -> None:
    """
    Generate reference database files based on configuration.
    
    Args:
        config_file (str): Path to the configuration file.
    """
    with CONSOLE.status("[bold green]Database generation started[/bold green]") as status:
        # Load and validate configuration
        config = load_config(config_file)
        status.update(
            status="[bold green]Configuration file parsed and validated successfully[/bold green]",
            spinner="squish",
            spinner_style="yellow"
        )
        
        if not validate_config(config):
            CONSOLE.log("Configuration validation failed", style="red")
            sys.exit(1)
            
        # Verify paths
        if not verify_paths(config):
            CONSOLE.log("Path verification failed", style="red")
            sys.exit(1)
            
        time.sleep(1)  # Visual feedback
        status.update(
            status="[bold green]Reference files checked and verified[/bold green]",
            spinner="squish",
            spinner_style="yellow"
        )
        
        # Process allotype data
        process_allotype_data(config)
        
        time.sleep(1)  # Visual feedback
        status.update(
            status="[bold green]HLA/MHC list generated and database output created successfully[/bold green]",
            spinner="squish",
            spinner_style="yellow"
        )
        
        time.sleep(1)  # Visual feedback
        status.update(
            status="[bold green]Database generation completed[/bold green]",
            spinner="squish",
            spinner_style="green"
        )