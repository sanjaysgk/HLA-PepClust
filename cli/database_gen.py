"""Verify that the database is generated correctly."""
import os
import sys
import json
import pandas as pd 
from cli.logger import CONSOLE
import time

def _prase_config_file(config_file):
    if not os.path.exists(config_file):
        # sys.exit(f"Config file {config_file} does not exist")
        CONSOLE.log(f"Config file {config_file} does not exist", style="red")
        CONSOLE.log("Using default configuration", style="blue")
        default_config ={
            "human":{
                "path": "Gibbs_motifs_human",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/hla_data.csv",
                "ref_data": "data/ref_data"
            },

            "mouse" :{
                "path": "Gibbs_motifs_mouse",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/mhc_data.csv",
                "ref_data": "data/ref_data"
            }
        }
        print(default_config)
        return default_config 
           
    else:
        with open(config_file, "r") as f:
            config = json.load(f)
        return config

def _check_ref_files(config):
    for species in config:
        if not os.path.exists(os.path.join(config[species]["ref_data"], config[species]["path"])):
            # sys.exit(f"Path {os.path.join(config[species]['ref_data'], config[species]['path'])} does not exist")
            # print("Path", config[species]["path"])
            # print("Ref data", config[species]["ref_data"])
            # print("Join", os.path.join(config[species]["ref_data"], config[species]["path"].lstrip('/')))
            CONSOLE.log(f"Path {os.path.join(config[species]['ref_data'],config[species]['path'])} does not exist. Check the path.", style="red")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]["ref_data"],config[species]["path"],config[species]["matrix"])):
            # sys.exit(f"Matrix path {config[species]['path'] + config[species]['matrix']} does not exist")
            CONSOLE.log(f"Matrix path {os.path.join(config[species]['ref_data'],config[species]['path'],config[species]['matrix'])} does not exist")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]["motif"])):
            # sys.exit(f"Motif path {config[species]['path'] + config[species]['motif']} does not exist")
            CONSOLE.log(f"Motif path {os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]['motif'])} does not exist")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]["ref_data"], config[species]["allotypes"])):
            # sys.exit(f"Allotype path {config[species]['path'] + config[species]['allotypes']} does not exist")
            CONSOLE.log(f"Allotype path {os.path.join(config[species]['ref_data'], config[species]['allotypes'])} does not exist")
            sys.exit(1)

def _HLA_liist(config):
    for species in config:
        allotype_file = os.path.join(config[species]["ref_data"], config[species]["allotypes"])
        allotypes = pd.read_csv(allotype_file)
        if species == "human":
            allotypes.rename(columns={ "formatted_HLA":"formatted_allotypes","HLA":"allotypes"}, inplace=True)
            allotypes['motif'] = "HLA_" + allotypes['formatted_allotypes'] + ".png"
            allotypes['species'] = species
        elif species == "mouse":
            allotypes.rename(columns={"formatted_MHC":"formatted_allotypes","MHC":"allotypes"}, inplace=True)
            allotypes['motif'] = allotypes['formatted_allotypes'] + ".png"
            allotypes['species'] = species
        else:
            # sys.exit(f"Species {species} not supported")
            CONSOLE.log(f"Species {species} not supported")
            sys.exit(1)
        allotypes['motif_path'] = ""
        allotypes['matrices_path'] = ""
        for motifs in allotypes['motif']:
            CONSOLE.log(f"[yellow]{species}[/yellow] {config[species]['path']}/motif/{motifs}", style="blue")
            if not os.path.exists(f"{config[species]['ref_data']}/{config[species]['path']}/motif/{motifs}"):
                # sys.exit(f"Motif {motifs} does not exist")
                CONSOLE.log(f"Motif {motifs} does not exist")
                sys.exit(1)
            elif not os.path.exists(f"{os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]['matrix'])}/{motifs.replace('.png', '.txt')}"):
                # sys.exit(f"Matrix {motifs.replace('.png', '.txt')} does not exist")
                CONSOLE.log(f"Matrix {motifs.replace('.png', '.txt')} does not exist")
                sys.exit(1)
            else:
                allotypes.loc[allotypes['motif'] == motifs, 'motif_path'] = f"{config[species]['path']}/motif/{motifs}"
                allotypes.loc[allotypes['motif'] == motifs, 'matrices_path'] = f"{config[species]['path']}/matrices/{motifs.replace('.png', '.txt')}"
                # CONSOLE.log(f"Motif {motifs} exists")
        # allotype_file_db = allotype_file.replace('.csv', '.db')        
        allotypes.to_csv(f"{config[species]['ref_data']}/{species}.db", index=False)
        CONSOLE.log(f"Database {species}.db created successfully ({config[species]['ref_data']}/{species}.db)", style="green")
        
def Database_gen(config_file):
    # CONSOLE.print("Database generation started", style="green")
    with CONSOLE.status("[bold green] Database generation started [bold green]") as status:
        config = _prase_config_file(config_file)
        # CONSOLE.log("Configuration file parsed and validated successfully.", style="green")
        status.update(
                        status=f"[bold green] Configuration file parsed and validated successfully. [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        
        _check_ref_files(config)
        time.sleep(3)
        
        # CONSOLE.log("Reference files checked and verified.", style="green")
        status.update(
                        status=f"[bold green] Reference files checked and verified [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        _HLA_liist(config)
        time.sleep(3)
        # CONSOLE.log("HLA list generated and database output created successfully.", style="green")
        status.update(
                        status=f"[bold green] HLA list generated and database output created successfully [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        time.sleep(3)
        # CONSOLE.print("Database generation completed", style="green")
        status.update(
                        status=f"[bold green] Database generation completed [/bold green]",
                        spinner="squish",
                        spinner_style="green",
                            )
        sys.exit(0)
    
# Database_gen("config.json")
# if __name__ == "__main__":
#     config_file = "config.json"
#     config = _prase_config_file(config_file)
#     _check_ref_files(config)
#     hla_list = _HLA_liist(config)
#     # print(hla_list)
#     # print(config)
#     CONSOLE.log("Config file parsed successfully")
#     sys.exit(0)