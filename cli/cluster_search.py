"""
Cluster Search module for HLA-PEPCLUST.

This module handles the core functionality of comparing peptide clusters with
reference HLA/MHC motifs to identify the best matches.
"""

import os
import sys
import time
import re
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import altair as alt
from typing import Dict, List, Tuple, Union, Optional
from PIL import Image, ImageDraw, ImageFont
from matplotlib.colors import LinearSegmentedColormap
from jinja2 import Template

from cli.logger import CONSOLE, display_search_results, save_console_log
from cli.html_config import html_content, body_start
from cli.utils import (
    generate_unique_random_ids, 
    ensure_directory_exists,
    format_hla_name,
    amino_acid_order_identical,
    format_input_gibbs,
    find_file_by_pattern,
    create_output_directories
)


class ClusterSearch:
    """
    Class for searching and correlating peptide clusters with reference HLA/MHC motifs.
    """
    
    def __init__(self):
        """Initialize the ClusterSearch class."""
        self.correlation_dict = {}
        self.valid_HLA = []
        self.console = CONSOLE
        self.species = None
        self.data_dir = None
        self._outfolder = None

    def _db_loader(self, db_path: str, species: str) -> pd.DataFrame:
        """
        Load the database file.

        Args:
            db_path (str): Path to the database file
            species (str): Species name ('human' or 'mouse')
            
        Returns:
            pd.DataFrame: DataFrame containing the database
            
        Raises:
            FileNotFoundError: If database file doesn't exist
        """
        species = str(species).lower()
        db_file = os.path.join(db_path, f'{species}.db')
        
        if not os.path.exists(db_file):
            raise FileNotFoundError(f"Database file {db_file} does not exist.")
        
        self.data_dir = db_path
        self.species = species
        
        return pd.read_csv(db_file)

    def _make_dir(self, path: str) -> str:
        """
        Create output directory for results.

        Args:
            path (str): Base path for the output directory
            
        Returns:
            str: Path to the created output directory
        """
        output_dir = os.path.join(path, 'clust_result')
        ensure_directory_exists(output_dir)
        
        self.console.log(
            f"Output directory created: {output_dir}", 
            style="blue"
        )
        
        return output_dir

    def check_HLA_DB(self, HLA_list: List[str], ref_folder: str) -> bool:
        """
        Check if the provided HLA types exist in the reference database.

        Args:
            HLA_list (List[str]): List of HLA types to check
            ref_folder (str): Path to the reference folder
            
        Returns:
            bool: True if all HLA types are valid, False otherwise
        """
        if not HLA_list:
            self.console.log(
                "No HLA types provided. Using all available HLA types from the reference folder."
            )
            return True

        # Get list of available HLA types from reference folder
        DB_hla_list = [
            format_hla_name(filename) for filename in os.listdir(ref_folder)
        ]

        # Check each provided HLA type
        for HLA in HLA_list:
            formatted_hla = format_hla_name(HLA)
            if formatted_hla not in DB_hla_list:
                self.console.print(
                    f"HLA type {HLA} not found in the reference folder."
                )
                return False
            else:
                self.valid_HLA.append(HLA)

        return True

    def compute_correlations_v2(
        self,
        db: pd.DataFrame,
        gibbs_results: str,
        n_clusters: str,
        output_path: str,
        hla_list: Optional[List[str]] = None,
    ) -> None:
        """
        Compute correlations between test and reference Gibbs matrices.

        Args:
            db (pd.DataFrame): Database DataFrame
            gibbs_results (str): Path to Gibbs results
            n_clusters (str): Number of clusters to use ('all', 'best_KL', or specific number)
            output_path (str): Path to save results
            hla_list (List[str], optional): List of HLA types to process
        """
        # Check for Gibbs matrices
        gibbs_result_matrix = os.path.join(gibbs_results, "matrices")
        if not (os.path.exists(gibbs_result_matrix) and 
                any(".mat" in file for file in os.listdir(gibbs_result_matrix))):
            raise FileNotFoundError(f"No Gibbs matrices found in {gibbs_result_matrix}")
            
        # Create output directory
        self._outfolder = self._make_dir(output_path)

        # Log HLA types being processed
        if hla_list:
            if not isinstance(hla_list, list):
                raise TypeError("HLA types must be provided as a list [Check main].")
            self.console.log(f"Processing specific HLA types: {hla_list}")
        else:
            hla_list = None

        # Start timing the process
        start_time = time.time()
        cluster_found = []
        
        # Process clusters based on requested number
        if n_clusters == "all":
            self._process_all_clusters(gibbs_result_matrix, db, hla_list)
        elif n_clusters == "best_KL":
            self._process_best_kl_clusters(gibbs_result_matrix, db, hla_list)
        elif n_clusters.isdigit() and 0 < int(n_clusters) <= 6:
            self._process_specific_clusters(gibbs_result_matrix, db, n_clusters, hla_list, cluster_found)
        else:
            self.console.log(
                f"Given n_clusters param {n_clusters} is invalid. Proceeding with 'all' clusters"
            )
            self._process_all_clusters(gibbs_result_matrix, db, hla_list)

        # Report processing time
        end_time = time.time()
        elapsed_time = end_time - start_time
        self.console.log(
            f"Cluster Search process completed in {elapsed_time:.2f} seconds."
        )
        
        # Check if any clusters were found (for specific cluster number)
        if n_clusters.isdigit() and 0 < int(n_clusters) <= 6 and len(cluster_found) == 0:
            self.console.log(
                f"No cluster files found for {n_clusters} clusters. Exiting.")
            sys.exit(1)

    def _process_all_clusters(self, gibbs_matrix_dir: str, db: pd.DataFrame, hla_list: Optional[List[str]]) -> None:
        """
        Process all clusters in the Gibbs matrix directory.
        
        Args:
            gibbs_matrix_dir (str): Directory containing Gibbs matrices
            db (pd.DataFrame): Database DataFrame
            hla_list (List[str], optional): List of HLA types to process
        """
        self.console.log("Processing all clusters", style="blue")
        with self.console.status("Processing all clusters") as status:
            for gibbs_f in os.listdir(gibbs_matrix_dir):
                for mat_path in db['matrices_path']:
                    # Check if we should process this matrix
                    should_process = False
                    if hla_list is not None:
                        hla_id = format_hla_name(str(mat_path).split('/')[0])
                        if hla_id in hla_list:
                            should_process = True
                    else:
                        should_process = True
                    
                    if should_process:
                        correlation = self._compute_and_log_correlation_V2(
                            os.path.join(gibbs_matrix_dir, gibbs_f),
                            mat_path,
                        )
                        status.update(
                            status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                            spinner="squish",
                            spinner_style="yellow",
                        )
                    else:
                        self.console.log(
                            f"Skipping {str(mat_path).split('/')[0]} as it is not in the provided HLA list"
                        )

    def _process_best_kl_clusters(self, gibbs_matrix_dir: str, db: pd.DataFrame, hla_list: Optional[List[str]]) -> None:
        """
        Process clusters with best KL divergence.
        
        Args:
            gibbs_matrix_dir (str): Directory containing Gibbs matrices
            db (pd.DataFrame): Database DataFrame
            hla_list (List[str], optional): List of HLA types to process
        """
        self.console.log("Processing for best KL divergence clusters")
        with self.console.status("Processing best KL divergence clusters") as status:
            for gibbs_f in os.listdir(gibbs_matrix_dir):
                for mat_path in db['matrices_path']:
                    # Check if we should process this matrix
                    should_process = False
                    if hla_list is not None:
                        hla_id = format_hla_name(str(mat_path).split('/')[0])
                        if hla_id in hla_list:
                            should_process = True
                    else:
                        should_process = True
                    
                    if should_process:
                        correlation = self._compute_and_log_correlation_V2(
                            os.path.join(gibbs_matrix_dir, gibbs_f),
                            mat_path,
                        )
                        status.update(
                            status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                            spinner="squish",
                            spinner_style="yellow",
                        )

    def _process_specific_clusters(self, gibbs_matrix_dir: str, db: pd.DataFrame, 
                                 n_clusters: str, hla_list: Optional[List[str]], 
                                 cluster_found: List[str]) -> None:
        """
        Process clusters for a specific number of clusters.
        
        Args:
            gibbs_matrix_dir (str): Directory containing Gibbs matrices
            db (pd.DataFrame): Database DataFrame
            n_clusters (str): Number of clusters to use
            hla_list (List[str], optional): List of HLA types to process
            cluster_found (List[str]): List to collect found cluster files
        """
        self.console.log(f"Processing for {n_clusters} clusters")
        with self.console.status(f"Processing {n_clusters} clusters") as status:
            for gibbs_f in os.listdir(gibbs_matrix_dir):
                if gibbs_f.endswith(f"of{n_clusters}.mat"):
                    cluster_found.append(gibbs_f)
                    for mat_path in db['matrices_path']:
                        # Check if we should process this matrix
                        should_process = False
                        if hla_list is not None:
                            hla_id = format_hla_name(str(mat_path).split('/')[0])
                            if hla_id in hla_list:
                                should_process = True
                        else:
                            should_process = True
                        
                        if should_process:
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_matrix_dir, gibbs_f),
                                mat_path,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )

    def _compute_and_log_correlation_V2(
        self,
        gibbs_f: str,
        ref_mat: str,
    ) -> float:
        """
        Compute correlation between a Gibbs matrix and a reference matrix.

        Args:
            gibbs_f (str): Path to the Gibbs matrix file
            ref_mat (str): Path to the reference matrix
            
        Returns:
            float: Correlation value or NaN if an error occurs
        """
        try:
            # Create output directory if it doesn't exist
            ensure_directory_exists(os.path.join(self._outfolder, 'corr-data'))
            
            # Format input data
            m1 = format_input_gibbs(gibbs_f)
            m2 = format_input_gibbs(os.path.join(self.data_dir, ref_mat))

            # Align amino acid order
            m1, m2 = amino_acid_order_identical(m1, m2)
            
            # Calculate correlation
            correlation = m1.corrwith(m2, axis=1).mean()
            
            # Store the result in the correlation dictionary
            self.correlation_dict[(gibbs_f), (ref_mat)] = correlation
            
            return correlation

        except Exception as e:
            self.console.print(
                f"Failed to compute correlation between {gibbs_f} and {ref_mat}: {str(e)}"
            )
            return float('nan')

    def find_highest_correlation_for_each_row(self, correlation_dict: Dict) -> Dict:
        """
        Find the highest correlation for each cluster matrix.

        Args:
            correlation_dict (Dict): Dictionary of correlation values
            
        Returns:
            Dict: Dictionary mapping each cluster to its best match and correlation
        """
        rows = sorted(set(key[0] for key in correlation_dict.keys()))
        cols = sorted(set(key[1] for key in correlation_dict.keys()))

        # Create a DataFrame to store correlation values
        matrix = pd.DataFrame(index=rows, columns=cols, dtype=float)

        # Populate the correlation matrix
        for (row, col), value in correlation_dict.items():
            matrix.loc[row, col] = value

        highest_corr_per_row = {}

        # Find column with highest correlation for each row
        for row in matrix.index:
            highest_col = matrix.loc[row].idxmax()  # Find column with highest correlation for the row
            highest_corr = matrix.loc[row, highest_col]  # Get the highest correlation value

            highest_corr_per_row[row] = (highest_col, highest_corr)

        return highest_corr_per_row

    def plot_heatmap(self, output_path: str) -> None:
        """
        Generate a heatmap visualization of correlation values.
        
        Args:
            output_path (str): Path to save the heatmap
        """
        # Check if correlation_dict is empty
        if not self.correlation_dict:
            raise ValueError(
                "correlation_dict is empty. Cannot generate heatmap. [Tip: try running without --n_clusters flag or check your gibbs output files]")

        rows = sorted(set(key[0] for key in self.correlation_dict.keys()))
        cols = sorted(set(key[1] for key in self.correlation_dict.keys()))

        # Create an empty DataFrame with the given rows and columns
        matrix = pd.DataFrame(index=rows, columns=cols, dtype=float)
        
        # Fill the matrix with the correlation values
        for (row, col), value in self.correlation_dict.items():
            matrix.loc[row, col] = value

        # Handle NaN values by filling them with 0
        matrix.fillna(0, inplace=True)

        # Check if the matrix is empty after filling NaN values
        if matrix.empty or matrix.isnull().all().all():
            raise ValueError(
                "Matrix is empty or filled with NaN values. Cannot generate heatmap. [Tip: try running without --n_clusters flag]")

        # Create custom colormap
        custom_cmap = LinearSegmentedColormap.from_list(
            "CustomColours",
            [
                "white", "white", "white", "white", "white", "white", "whitesmoke",
                "lightgrey", "grey", "deeppink",
            ],
        )

        # Count HLA types for dividing lines
        HLA_A_count = sum(col.startswith("HLA_A") for col in matrix.columns)
        HLA_B_count = sum(col.startswith("HLA_B") for col in matrix.columns)

        # Create the figure and axis
        fig, ax = plt.subplots(figsize=(20, 7.5))

        # Plot the heatmap
        sns.heatmap(
            matrix,
            annot=False,
            cmap=custom_cmap,
            cbar=True,
            linewidths=0.5,
            square=True,
            ax=ax,
        )

        # Add vertical and horizontal lines
        y_min, y_max = ax.get_ylim()
        x_min, x_max = ax.get_xlim()

        plt.title("Correlation Heatmap")
        plt.xlabel("HLA Reference")
        plt.ylabel("Input Samples")
        plt.xticks(rotation=90, ha="center")

        # Add division lines for different HLA types
        ax.vlines(HLA_A_count, y_min, y_max, color="k")
        ax.vlines(HLA_A_count + HLA_B_count, y_min, y_max, color="k")

        # Add division lines for different sample clusters
        unique_names = matrix.index.unique()
        base_names = pd.Series([re.sub(r"_gibbs\..*", "", str(name)) for name in unique_names])
        occurrence_counts = {base: sum(base_names.str.startswith(base)) for base in base_names.unique()}

        starting_horizontal = 0
        for count in occurrence_counts.values():
            starting_horizontal += count
            ax.hlines(
                starting_horizontal,
                x_min,
                x_max,
                color="grey",
                linewidth=1,
                linestyles="--",
            )

        # Save the heatmap
        plt.tight_layout()
        plt.savefig(f"{self._outfolder}/heatmap.png")
        plt.close()

    def _make_correlation_plot(self, gibbs_df: pd.DataFrame, motif_df: pd.DataFrame, mat_motif: str) -> None:
        """
        Create a correlation plot comparing amino acid frequencies.
        
        Args:
            gibbs_df (pd.DataFrame): Gibbs matrix DataFrame
            motif_df (pd.DataFrame): Reference motif DataFrame
            mat_motif (str): Identifier for the output file
        """
        # Prepare DataFrames
        df1 = pd.DataFrame(gibbs_df)
        df2 = pd.DataFrame(motif_df)
        
        df1['Position'] = df1.index + 1
        df2['Position'] = df2.index + 1
        
        df1_melted = df1.melt(id_vars='Position', var_name='Amino Acid', value_name='PWM1')
        df2_melted = df2.melt(id_vars='Position', var_name='Amino Acid', value_name='PWM2')

        # Merge the two dataframes
        df_merged = pd.merge(df1_melted, df2_melted, on=['Position', 'Amino Acid'])
        
        # Create the base chart
        base = alt.Chart(df_merged, width="container").mark_circle().encode(
            x='PWM1',
            y='PWM2',
            color='Amino Acid',
            tooltip=['Amino Acid', 'PWM1', 'PWM2', 'Position']
        )
        
        # Calculate the correlation coefficient
        corr_coef = df_merged[['PWM1', 'PWM2']].corr().iloc[0, 1]

        # Create the regression line
        regression_line = base.transform_regression('PWM1', 'PWM2').mark_line(opacity=0.50, shape='mark').transform_fold(
            ["reg-line"], 
            as_=["Regression", "y"]
        ).encode(alt.Color("Regression:N"))
        
        # Add the correlation coefficient as text
        corr_text = alt.Chart(pd.DataFrame({
            'PWM1': [df_merged['PWM1'].min()],
            'PWM2': [df_merged['PWM2'].max()],
            'text': [f'Correlation: {corr_coef:.2f}']
        })).mark_text(align='left', baseline='top', dx=8, dy=-5).encode(
            x='PWM1:Q',
            y='PWM2:Q',
            text='text:N'
        )
        
        # Combine the charts
        chart = base + regression_line + corr_text
        
        # Ensure the output directory exists
        ensure_directory_exists(os.path.join(self._outfolder, 'corr-data'))
        
        # Save as JSON and PNG
        chart.save(f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.json")
        chart.save(f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.png")

    def insert_script_png_json(self, script_data_path: str, img_fallback_path: str, div_id: str) -> str:
        """
        Generate JavaScript for loading a Vega visualization with PNG fallback.
        
        Args:
            script_data_path (str): Path to the JSON data file
            img_fallback_path (str): Path to the fallback image
            div_id (str): HTML div ID for the visualization
            
        Returns:
            str: JavaScript code
        """
        script_template = Template('''
            fetch('{{ script_data_path }}')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", '{{ script_data_path }}', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#{{ div_id }}", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("{{ div_id }}").innerHTML = `<img src="{{ img_fallback_path }}" alt="Fallback Image" width="100%">`;
                });
        ''')

        return script_template.render(
            script_data_path=script_data_path, 
            img_fallback_path=img_fallback_path, 
            div_id=div_id
        )

    def render_hla_section(self, hla_name: str, correlation_chart_id: str, 
                          best_cluster_img: str, naturally_presented_img: str) -> str:
        """
        Generate HTML for an HLA section.
        
        Args:
            hla_name (str): Name of the HLA
            correlation_chart_id (str): ID for the correlation chart
            best_cluster_img (str): Path to the best cluster image
            naturally_presented_img (str): Path to the naturally presented image
            
        Returns:
            str: HTML code for the section
        """
        template = Template('''
        <div class="row" style="border: 2px solid #007bff;">
        <div class="row">
            <h3 style="text-align: center;">{{ hla_name }}</h3>
        </div>
        <div class="row">
            <div class="col">
            <div class="card">
                <h5 style="text-align: center;">Correlation</h5>
                <div id="{{correlation_chart_id}}"></div>
            </div>
            </div>
            <div class="col">
            <div class="card">
                <h5 style="text-align: center;">Identified Best Cluster</h5>
                <div class="card">
                {% if best_cluster_img and best_cluster_img != "None" %}
                <img src="{{ best_cluster_img }}" class="bd-placeholder-img card-img" width="100%" height="260" alt="No Image">
                {% else %}
                <div class="no-img" width="100%" height="260"></div>
                {% endif %}
                <div style="position: absolute; top: 5px; right: 10px;">
                    <i class="bi bi-three-dots-vertical" id="dropdownMenuButton" data-bs-toggle="dropdown" aria-expanded="false"></i>
                    <div class="dropdown-menu">
                    <h6 class="dropdown-header">Download</h6>
                    {% if best_cluster_img and best_cluster_img != "None" %}
                    <a class="dropdown-item" href="{{ best_cluster_img }}">png</a>
                    {% endif %}
                    </div>
                </div>
                </div>
            </div>
            </div>
            <div class="col">
            <div class="card">
                <h5 style="text-align: center;">Naturally presented {{ hla_name }}</h5>
                <div class="card">
                {% if naturally_presented_img and naturally_presented_img != "None" %}
                <img src="{{ naturally_presented_img }}" class="bd-placeholder-img card-img" width="100%" height="260" alt="No Image">
                {% else %}
                <div class="no-img"></div>
                {% endif %}
                <div style="position: absolute; top: 5px; right: 10px;">
                    <i class="bi bi-three-dots-vertical" id="dropdownMenuButton" data-bs-toggle="dropdown" aria-expanded="false"></i>
                    <div class="dropdown-menu">
                    <h6 class="dropdown-header">Download</h6>
                    {% if naturally_presented_img and naturally_presented_img != "None" %}
                    <a class="dropdown-item" href="{{ naturally_presented_img }}">Download</a>
                    {% endif %}
                    </div>
                </div>
                </div>
            </div>
            </div>
        </div>
        </div>
        ''')
        
        return template.render(
            hla_name=hla_name, 
            correlation_chart_id=correlation_chart_id, 
            best_cluster_img=best_cluster_img, 
            naturally_presented_img=naturally_presented_img
        )

    def make_datatable(self, correlation_dict: Dict) -> pd.DataFrame:
        """
        Convert correlation dictionary to a DataFrame.
        
        Args:
            correlation_dict (Dict): Dictionary of correlation values
            
        Returns:
            pd.DataFrame: Formatted DataFrame
        """
        df = pd.DataFrame(correlation_dict.items(), columns=['HLA', 'Correlation'])
        df['Cluster'] = df['HLA'].apply(lambda x: x[0].split('/')[-1].split('.')[1])    
        df['HLA'] = df['HLA'].apply(lambda x: x[1])
        df['HLA'] = df['HLA'].apply(lambda x: x.split('/')[-1].replace('.txt', ''))
        df['Correlation'] = df['Correlation'].apply(lambda x: round(x, 2))
        df = df.sort_values(by='Correlation', ascending=False)
        df = df.reset_index(drop=True)
        df = df[['Cluster', 'HLA', 'Correlation']]
        return df

    def make_heatmap(self, correlation_dict: Dict, threshold: float = 0.5) -> bool:
        """
        Generate an interactive heatmap visualization.
        
        Args:
            correlation_dict (Dict): Dictionary of correlation values
            threshold (float): Threshold for color scaling
            
        Returns:
            bool: True if successful, False otherwise
        """
        df = self.make_datatable(correlation_dict)
        
        try:
            # Save CSV
            df.to_csv(os.path.join(self._outfolder, 'corr-data', 'corr_matrix.csv'), index=False)
            
            # Create heatmap chart
            chart_h = alt.Chart(df, width="container").mark_rect().encode(
                alt.X("HLA:O").title("HLA").axis(labelAngle=-45),
                alt.Y("Cluster:O").title("Cluster"),
                alt.Color(
                    "Correlation",
                    scale=alt.Scale(
                        domain=[df["Correlation"].min(), threshold, df["Correlation"].max()],
                        range=["#d3d3d3", "#add8e6", "#ff4500"],  # Low values fade, high values bright
                    ),
                    legend=None
                ),
                tooltip=[
                    alt.Tooltip("HLA", title="HLA"),
                    alt.Tooltip("Cluster", title="Cluster"),
                    alt.Tooltip("Correlation", title="Correlation"),
                ],
            ).configure_view(
                step=13,
                strokeWidth=0
            ).configure_axis(
                domain=False
            )
            
            # Save the chart
            chart_h.save(f"{os.path.join(self._outfolder,'corr-data')}/correlation_heatmap.json")
            chart_h.save(f"{os.path.join(self._outfolder,'corr-data')}/correlation_heatmap.png")
            return True
            
        except Exception as e:
            self.console.log(f"Failed to save the correlation heatmap: {str(e)}")
            return False

    def make_datatable_html(self, correlation_dict: Dict, df: Optional[pd.DataFrame] = None, threshold: float = 0.5) -> str:
        """
        Generate HTML table for correlation data.
        
        Args:
            correlation_dict (Dict): Dictionary of correlation values
            df (pd.DataFrame, optional): Pre-processed DataFrame
            threshold (float): Threshold for highlighting high correlations
            
        Returns:
            str: HTML table code
        """
        if df is None:
            df = self.make_datatable(correlation_dict)
            
        table_start = """
        <table id="correlation_table" class="table table-bordered">
          <thead class="thead-dark">
            <tr>
              <th scope="col">ID's</th>
              <th scope="col">Cluster</th>
              <th scope="col">Best HLA/MHC Match</th>
              <th scope="col">Score</th>
              <th scope="col">KLD Score</th>
            </tr>
          </thead>
          <tbody>
        """
        
        tr = ""
        for i in range(len(df)):
            if df['Correlation'][i] >= threshold:
                tr += f"""
                <tr class="table-success">
                    <td>{i+1}</td>
                    <td>{df['Cluster'][i]}</td>
                    <td>{df['HLA'][i]}</td>
                    <td>{df['Correlation'][i]}</td>
                    <td>0.0</td>
                </tr>
                """
            elif df['Correlation'][i] >= 0.5 and df['Correlation'][i] < threshold:
                tr += f"""
                <tr class="table-warning">
                    <td>{i+1}</td>
                    <td>{df['Cluster'][i]}</td>
                    <td>{df['HLA'][i]}</td>
                    <td>{df['Correlation'][i]}</td>
                    <td>0.0</td>
                </tr>
                """
            else:
                tr += f"""
                <tr class="table-danger">
                    <td>{i+1}</td>
                    <td>{df['Cluster'][i]}</td>
                    <td>{df['HLA'][i]}</td>
                    <td>{df['Correlation'][i]}</td>
                    <td>0.0</td>
                </tr>
                """
                
        table_end = """
           </tbody>
        </table>
        """
        
        return table_start + tr + table_end

    def generate_html_layout(self, correlation_dict: Dict, db: pd.DataFrame, gibbs_out: str, immunolyser: bool = False) -> None:
        """
        Generate HTML report layout with all visualization components.
        
        Args:
            correlation_dict (Dict): Dictionary of correlation values
            db (pd.DataFrame): Database DataFrame
            gibbs_out (str): Path to Gibbs output
            immunolyser (bool): Whether to generate immunolyser output
        """
        # Find highest correlation for each cluster
        highest_corr_per_row = self.find_highest_correlation_for_each_row(correlation_dict)
        display_search_results(highest_corr_per_row, 0.8)

        # Create output directories
        dirs = create_output_directories(self._outfolder)
        
        output_dict = {}
        html_create = html_content + body_start
        
        # JavaScript for end of body
        body_end_1 = """
</div>
</div>
</div>
</main>
  <!-- Datatable  -->
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.2/dist/js/bootstrap.bundle.min.js"
    integrity="sha384-MrcW6ZMFYlzcLA8Nl+NtUVF0sA7MsXsP1UyJoMp4YLEuNSfAP+JcXn/tWtIaxVXM"
    crossorigin="anonymous"></script>
  <script src="https://code.jquery.com/jquery-3.7.1.js"></script>
  <script src="https://cdn.datatables.net/2.2.2/js/dataTables.js"></script>
  <script src="https://cdn.datatables.net/2.2.2/js/dataTables.bootstrap5.js"></script>

<footer class="text-body-secondary py-5">
    <div class="container">
      <p class="float-end mb-1">
        <a href="#">Back to top</a>
      </p>
      <p class="mb-1"> Purcell Lab, Monash University&copy; 2025, please refer to git repo</p>
      <p class="mb-0">Please visit our <a href="https://github.com/Sanpme66/HLA-PepClust">GitHub repository</a> or read our <a
          href="https://github.com/Sanpme66/HLA-PepClust">getting started guide</a>.</p>
    </div>
  </footer>
  <script>
   $(document).ready(function() {
    $('#correlation_table').DataTable();
  });
  document.addEventListener('DOMContentLoaded', function() {
    const downloadLinks = document.querySelectorAll('.dropdown-item');

    downloadLinks.forEach(link => {
      link.addEventListener('click', function(event) {
        event.preventDefault();
        const url = this.getAttribute('href');
        const a = document.createElement('a');
        a.href = url;
        a.download = url.split('/').pop();
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
      });
    });
  });
  document.addEventListener("DOMContentLoaded", function () {
    const noImgDivs = document.querySelectorAll(".no-img"); // Select all elements with class "no-img"
    const placeholderSVG = `
        <svg class="bd-placeholder-img card-img" width="100%" height="260" xmlns="http://www.w3.org/2000/svg"
            role="img" aria-label="Placeholder: No Image" preserveAspectRatio="xMidYMid slice" focusable="false">
            <title>Placeholder</title>
            <rect width="100%" height="100%" fill="#868e96" />
            <text x="50%" y="50%" fill="#dee2e6" dy=".3em" text-anchor="middle">No image found \n check your gibbs output</text>
        </svg>
    `;
    noImgDivs.forEach(div => {
        div.innerHTML = placeholderSVG; // Add placeholder to each "no-img" div
    });
});
</script>

<script>
"""

        # Line break tag
        br_tag = """
        <br>
"""

        # End of HTML document
        body_end_2 = """
</script>
</body>
</html>
"""

        # JavaScript for heatmap
        heatmap_js = """
                // heatmap json or PNG loading
                fetch('corr-data/correlation_heatmap.json')
                .then(response => {
                    if (!response.ok) {
                        throw new Error("JSON fetch failed");
                    }
                    return response.json();
                })
                .then(data => {
                    console.log("Loaded JSON:", 'corr-data/correlation_heatmap.json', data);
                    var opt = {"renderer": "canvas", "actions": false};
                    vegaEmbed("#correlation_heatmap", data, opt);
                })
                .catch(error => {
                    console.error("Error fetching JSON, loading fallback image instead:", error);
                    document.getElementById("correlation_heatmap").innerHTML = `<img src="corr-data/correlation_heatmap.png" alt="Fallback Image" width="100%">`;
                });
        """

        # HTML for heatmap div
        heatmap_div = """
            <div class="row">
            <div class="col">
                <div class="card">
                    <h5 style="text-align: center;">Correlation Matrix</h5>
                    <div class="row">
                    <div id="correlation_heatmap"></div>
                    </div>
                </div>
            </div>
            </div>
            """

        # Generate heatmap JSON
        heatmap_json = self.make_heatmap(correlation_dict)
        
        # Generate HTML table
        df_corr_html = self.make_datatable_html(correlation_dict)
        
        # Prepare immunolyser output 
        immunolyser_out = ""
        immunolyser_out_js = ""
        
        # Process each cluster and its best match
        for row, (col, corr) in highest_corr_per_row.items():
            # Extract cluster ID from full path
            cluster_id = str(row).split("/")[-1].split('.')[1]
            
            # Initialize dictionary entry for this cluster
            output_dict[cluster_id] = {
                'cluster': str(row).split("/")[-1],
                'HLA': str(col).split('/')[-1].replace('.txt', ''),
                'correlation': corr,
                'gibbs_img': None,
                'nat_img': None,
                'corr_plot': None,
                'corr_json': None
            }
            
            # Find and copy cluster image
            gibbs_img = find_file_by_pattern(os.path.join(gibbs_out, 'logos'), cluster_id)
            if gibbs_img:
                dest_path = os.path.join(self._outfolder, 'cluster-img', f"{os.path.basename(gibbs_img)}")
                shutil.copy(gibbs_img, dest_path)
                output_dict[cluster_id]['gibbs_img'] = dest_path
            
            # Extract HLA ID based on species
            if self.species == 'human':
                hla = str(col).split('/')[-1].replace('.txt', '').split('_')[1]
            else:
                hla = str(col).split('/')[-1].replace('.txt', '')
                
            # Find and copy motif image
            nat_path = db[db['formatted_allotypes'] == hla]['motif_path'].values[0]
            if nat_path:
                dest_path = os.path.join(self._outfolder, 'allotypes-img', f"{os.path.basename(nat_path)}")
                shutil.copy(nat_path, dest_path)
                output_dict[cluster_id]['nat_img'] = dest_path
            
            # Generate correlation plot
            try:
                gibbs_mt = format_input_gibbs(row)
                nat_mat = format_input_gibbs(col)

                # Align amino acid order
                gibbs_mt, nat_mat = amino_acid_order_identical(gibbs_mt, nat_mat)
                
                # Create a unique identifier for this comparison
                mat_motif = f"{cluster_id}_{str(col).split('/')[-1].replace('.txt', '')}"
                
                # Generate the correlation plot
                self._make_correlation_plot(gibbs_mt, nat_mat, mat_motif)
                
                # Store paths to plot files
                output_dict[cluster_id]['corr_plot'] = f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.png"
                output_dict[cluster_id]['corr_json'] = f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.json"
            except Exception as e:
                self.console.print(
                    f"Failed to generate correlation plot for {row} and {col}: {str(e)}"
                )
            
            # Generate HTML for this HLA section
            rows_list = self.render_hla_section(
                hla, 
                f"correlation_chart_{cluster_id}", 
                str(output_dict[cluster_id]['gibbs_img']).replace(f"{self._outfolder}/", ''), 
                str(output_dict[cluster_id]['nat_img']).replace(f"{self._outfolder}/", '')
            )
            
            # Add to main HTML
            html_create += rows_list
            
            # Add JavaScript for the correlation chart
            plot_js = self.insert_script_png_json(
                str(output_dict[cluster_id]['corr_json']).replace(f"{self._outfolder}/", ''),
                str(output_dict[cluster_id]['corr_plot']).replace(f"{self._outfolder}/", ''),
                f"correlation_chart_{cluster_id}"
            )
            body_end_1 += plot_js
            
            # If immunolyser output is requested, add to those variables too
            immunolyser_out += rows_list
            immunolyser_out_js += plot_js
        
        # Add heatmap if available
        if heatmap_json:
            body_end_1 += heatmap_js
            immunolyser_out_js += heatmap_js
        else:
            heatmap_div = ""
            heatmap_js = ""
        
        # Complete the HTML
        html_create += br_tag + df_corr_html + br_tag + heatmap_div + body_end_1 + body_end_2
        
        # If immunolyser output is requested, generate those files
        if immunolyser:
            immunolyser_out += br_tag + br_tag + df_corr_html + heatmap_div
            with open(os.path.join(self._outfolder, "immunolyser-out.html"), "w") as file:
                file.write(immunolyser_out)
            with open(os.path.join(self._outfolder, "immunolyser-out.js"), "w") as file:
                file.write(immunolyser_out_js)
        
        # Save the main HTML file
        with open(os.path.join(self._outfolder, "clust-search-result.html"), "w") as file:
            file.write(html_create)
            
        self.console.log(f"HTML report saved to {os.path.join(self._outfolder, 'clust-search-result.html')}", style="bold green")


def run_cluster_search(args):
    """
    Main function to run the cluster search pipeline.
    
    Args:
        args: Command-line arguments from argparse
    """
    # Initialize ClusterSearch
    cluster_search = ClusterSearch()
    
    # Stage 1: Data processing
    cluster_search.console.rule(
        "[bold red]Stage 1/2: Data processing for correlation matrices."
    )
    
    # Validate species
    if str(args.species).lower() in ["mouse", "human"]:
        CONSOLE.log(
            f"Species provided: [bold yellow]{args.species}", style="bold green")
        CONSOLE.log(
            f"Loading reference database for [bold yellow]{args.species}", style="bold green")
        
        # Load database
        db = cluster_search._db_loader(args.reference_folder, args.species)
        CONSOLE.log(f"Reference database loaded successfully.",
                    style="bold green")
    else:
        raise ValueError(
            "Invalid species provided. Please provide a valid species. [Human or Mouse]")

    # Process HLA types if provided
    if args.hla_types:
        CONSOLE.log(
            f"HLA/MHC allotypes provided: [bold yellow]{args.hla_types}", style="bold green")
        u_hla_list = []
        
        with CONSOLE.status("Checking HLA/MHC in database") as status:
            status.update(
                status=f"[bold blue] Loading HLA/MHC in database",
                spinner="squish",
                spinner_style="yellow",
            )
            
            # If args.hla_types is a string, split it
            if isinstance(args.hla_types, str):
                hla_types = args.hla_types.split(',')
            else:
                # If it's already a list or other iterable, join and split to ensure proper format
                hla_types = ','.join(args.hla_types).split(',')
                
            # Check each HLA type against the database
            for u_hla in hla_types:
                if u_hla in db['formatted_allotypes'].values:
                    status.update(
                        status=f"[bold blue] HLA/MHC {u_hla} found in database",
                        spinner="squish",
                        spinner_style="yellow",
                    )
                    u_hla_list.append(u_hla)
                    CONSOLE.log(
                        f"HLA/MHC [bold yellow]{u_hla}[yellow] found in database", style="bold green")
                else:
                    status.update(
                        status=f"[bold blue] HLA/MHC {u_hla} not found in database",
                        spinner="squish",
                        spinner_style="yellow",
                    )
                    CONSOLE.log(
                        f"HLA/MHC [bold yellow]{u_hla} not found in database", style="bold red")
        
        # Use validated HLA list
        if len(u_hla_list) > 0:
            CONSOLE.log(
                f"Valid HLA/MHC allotypes provided: [bold yellow]{u_hla_list}", style="bold green")
            args.hla_types = u_hla_list
        else:
            args.hla_types = None
            CONSOLE.log(
                f"No valid HLA/MHC allotypes found. Using all available types.", style="bold yellow")
    else:
        args.hla_types = None
        CONSOLE.log(
            f"No HLA/MHC allotypes provided. Using all available types from the reference folder.", style="bold yellow")

    # Compute correlations
    CONSOLE.log("Calculating correlations between matrices...", style="bold blue")
    cluster_search.compute_correlations_v2(
        db,
        args.gibbs_folder,
        args.n_clusters,
        args.output,
        args.hla_types,
    )

    # Stage 2: Result visualization
    cluster_search.console.rule(
        "[bold red]Stage 2/2: Finding best matching Naturally presented HLA."
    )
    
    # Generate heatmap
    cluster_search.plot_heatmap(args.output)
    
    # Generate HTML report
    cluster_search.generate_html_layout(
        cluster_search.correlation_dict, 
        db, 
        args.gibbs_folder,
        args.immunolyser
    )
    
    # Save log if requested
    if args.log:
        log_file_path = os.path.join(cluster_search._outfolder, "search_cluster.log")
        save_console_log(log_file_path)