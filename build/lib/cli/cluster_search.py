from cli.imagegrid import imagelayout
from cli.logger import *
from cli.html_config import *
from cli.utils import *
from cli import __version__
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex
import logging
import argparse
from multiprocessing import Pool
from PIL import Image, ImageDraw, ImageFont
import os
import re
import time
from pathlib import Path
import json
from rich.traceback import install
from jinja2 import Template
import sys
import altair as alt

# import HTML
import shutil


install(show_locals=True)


# Configure logging
# logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# logging = CONSOLE
# _logConfig(logSave=True)
class ClusterSearch:
    def __init__(self):
        self.correlation_dict = {}
        self.valid_HLA = []
        self.console = CONSOLE
        self.species = None
        self.uhla = None
        self.data_dir = None
        # pd.DataFrame(columns=['Cluster', 'HLA', 'Correlation'])
        self.corr_df = None
        self.threshold = 0.70
        self.d3js_json = None
        self.db = None  # databse of motif and matrix
        self.treshold_img = False
        self.gibbs_results = None
        self.kld_df = None
   
        
    def generate_unique_random_ids(self, count: int) -> list:
        """
        Generate a list of unique 6-digit random IDs.

        :param count: Number of unique random IDs to generate
        :return: List of unique 6-digit random IDs
        """
        start = 100000
        end = 999999
        if count > (end - start + 1):
            raise ValueError(
                "Count is larger than the range of unique IDs available.")

        return np.random.choice(
            range(start, end + 1), size=count, replace=False
        ).tolist()

    def _db_loader(self, db_path: str, species: str) -> pd.DataFrame:
        """
        Load the database file.

        :param db_path: Path to the database file
        :return: DataFrame containing the database
        """
        self.db_path = db_path
        if not os.path.exists(os.path.join(db_path, f'{str(species).lower()}.db')):
            raise FileNotFoundError(f"Database file {db_path} does not exist.")
        else:
            self.data_dir = db_path
            self.species = str(species).lower()

        return pd.read_csv(os.path.join(db_path, f'{species}.db'))

    def parse_gibbs_output(self, gibbs_path: str, n_clusters: int) -> pd.DataFrame:
        """
        Parse the Gibbs output files.

        :param files_path: Path to the Gibbs output files
        :param n_clusters: Number of clusters
        :return: DataFrame with parsed data
        """
        res_path = os.path.join(gibbs_path, "res")
        if not os.path.exists(res_path):
            raise FileNotFoundError(f"Directory {res_path} does not exist.")

        for c_file in os.listdir(res_path):
            if f'{n_clusters}g.ds' in c_file:
                file_path = os.path.join(res_path, c_file)
                df = pd.read_csv(file_path, sep='\s+')
                # output_path = f'data/sampledata_701014/res_{n_clusters}g.csv'
                # df.to_csv(output_path, index=False)
                logging.info(f"Data parsed for No clusters {n_clusters}")
                return df

        raise FileNotFoundError(
            f"No cluster file found for {n_clusters} clusters.")

    def _make_dir(self, path: str, rand_ids: int) -> None:
        """
        Make a directory if it does not exist.

        :param path: Path to the directory
        """
        # if not os.path.exists(os.path.join(path, f'clust_result_{rand_ids}')):
        #     os.makedirs(os.path.join(path, f'clust_result_{rand_ids}'))

        if not os.path.exists(os.path.join(path, f'clust_result')):
            os.makedirs(os.path.join(path, f'clust_result'))
        # return os.path.join(path, f'clust_result_{rand_ids}')

        self.console.log(
            f"Output directory created {os.path.join(path, f'clust_result')}", style="blue")

        return os.path.join(path, f'clust_result')

    def _check_HLA_DB(self, HLA_list: list, ref_folder: str) -> bool:
        """
        Check if the HLA type is in the list.

        :param HLA_list: List of HLA types
        :param ref_folder: Reference folder containing HLA Matrix files
        :return: True if HLA type is in the list, False otherwise
        """
        if not HLA_list:
            # logging.warning("No HLA types provided. Using all available HLA types from the reference folder.")
            self.console.log(
                "No HLA types provided. Using all available HLA types from the reference folder."
            )
            return True

        DB_hla_list = [
            self.formate_HLA_DB(filename) for filename in os.listdir(ref_folder)
        ]
        # breakpoint()

        for HLA in self.formate_HLA_user_in(HLA_list):
            if self.formate_HLA_DB(HLA) not in DB_hla_list:
                # logging.error(f"HLA type {HLA} not found in the reference folder.")
                self.console.log(
                    f"HLA type {HLA} not found in the reference folder.")
                # logging.error(f"Available HLA types: {DB_hla_list}")
                self.console.log(f"Available HLA types: {DB_hla_list}")
                return False
            else:
                self.valid_HLA.append(HLA)

        return True

    @staticmethod
    def format_input_gibbs(gibbs_matrix: str) -> pd.DataFrame:
        """
        Format the Gibbs output matrix.

        :param gibbs_matrix: Path to the Gibbs matrix file
        :return: Processed DataFrame
        """
        df = pd.read_csv(gibbs_matrix)
        amino_acids = df.iloc[0, 0].split()
        df.iloc[:, 0] = df.iloc[:, 0].str.replace(
            r"^\d+\s\w\s", "", regex=True)
        new_df = pd.DataFrame(df.iloc[1:, 0].str.split(
            expand=True).values, columns=amino_acids)
        new_df.reset_index(drop=True, inplace=True)
        new_df = new_df.apply(pd.to_numeric, errors='coerce')
        return new_df

    @staticmethod
    def amino_acid_order_identical(df1: pd.DataFrame, df2: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        if list(df1.columns) != list(df2.columns):
            logging.warning(
                "The amino acid column order is different. Reordering columns.")
            df2 = df2[df1.columns]
        return df1, df2

    @staticmethod
    def formate_HLA_DB(HLA: str) -> str:
        """
        Format the HLA database file name to extract meaningful parts.

        :param HLA: HLA file name
        :return: Formatted HLA type
        """
        return HLA.replace("HLA_", "").replace("*", "").replace("txt", "").replace(".", "")

    def check_HLA_DB(self, HLA_list: list, ref_folder: str) -> bool:
        """
        Check if the HLA type is in the list.

        :param HLA_list: List of HLA types
        :param ref_folder: Reference folder containing HLA files
        :return: True if HLA type is in the list, False otherwise
        """
        if not HLA_list:
            # logging.warning("No HLA types provided. Using all available HLA types from the reference folder.")
            self.console.log(
                "No HLA types provided. Using all available HLA types from the reference folder."
            )
            return True

        DB_hla_list = [
            self.formate_HLA_DB(filename) for filename in os.listdir(ref_folder)
        ]

        for HLA in HLA_list:
            if self.formate_HLA_DB(HLA) not in DB_hla_list:
                # logging.error(f"HLA type {HLA} not found in the reference folder.")
                self.console.print(
                    f"HLA type {HLA} not found in the reference folder.")
                # logging.error(f"Available HLA types: {DB_hla_list}")
                # self.console.print(f"Available HLA types: {DB_hla_list}")
                return False
            else:
                self.valid_HLA.append(HLA)

        return True

    def compute_correlations(
        self,
        db: pd.DataFrame,
        gibbs_folder: str,
        human_reference_folder: str,
        n_clusters: str,
        output_path: str,
        hla_list: str = None,
    ) -> None:
        """
        Compute correlations between test and reference Gibbs matrices.

        :param gibbs_folder: Path to test matrices
        :param human_reference_folder: Path to reference matrices
        :param hla_list: List of HLA types to process (optional, default is all available types)
        """
        gibbs_matrix_folder = os.path.join(gibbs_folder, "matrices")
        if os.path.exists(gibbs_matrix_folder) and any(".mat" in file for file in os.listdir(gibbs_matrix_folder)):
            self.console.log(
                f"Found Gibbs matrices in {gibbs_matrix_folder}")
        else:
            raise FileNotFoundError(
                f"No Gibbs matrices found in {gibbs_matrix_folder}")
        # print(db.head())
        self._outfolder = self._make_dir(
            output_path, self.generate_unique_random_ids(6)[0]
        )  # self.generate_unique_random_ids(6)[0]

        if hla_list:
            hla_list = hla_list[0].split(
                ","
            )  # Split if comma-separated string is passed
            assert isinstance(
                hla_list, list), "HLA types must be provided as a list."
            self.console.log(f"Processing specific HLA types: {hla_list}")
        else:
            hla_list = None
            self.console.log("Processing all available HLA types.")

        # Check for valid HLA types in the reference folder
        if not self.check_HLA_DB(hla_list, human_reference_folder):
            self.console.log(
                "Invalid or missing HLA types. Aborting correlation computation."
            )
            return None
        # breakpoint()

        start_time = time.time()
        cluster_found = []

        if n_clusters == "all":
            self.console.log("Processing all clusters")
            with self.console.status("Processing all clusters") as status:
                for filename1 in os.listdir(gibbs_matrix_folder):
                    for filename2 in os.listdir(human_reference_folder):
                        # print(filename1, filename2,"####"*100)
                        if (
                            hla_list is None
                            or self.formate_HLA_DB(filename2) in hla_list
                        ):
                            correlation = self._compute_and_log_correlation(
                                gibbs_matrix_folder,
                                human_reference_folder,
                                filename1,
                                filename2,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {filename1} and {filename2} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )
                            # print(correlation)

        elif n_clusters == "best_KL":
            self.console.log("Processing for best KL divergence clusters")
            with self.console.status(
                "Processing best KL divergence clusters"
            ) as status:
                for filename1 in os.listdir(gibbs_matrix_folder):
                    for filename2 in os.listdir(human_reference_folder):
                        if (
                            hla_list is None
                            or self.formate_HLA_DB(filename2) in hla_list
                        ):
                            correlation = self._compute_and_log_correlation(
                                gibbs_matrix_folder,
                                human_reference_folder,
                                filename1,
                                filename2,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {filename1} and {filename2} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )

        elif n_clusters.isdigit() and 0 < int(n_clusters) <= 6:
            self.console.log(f"Processing for {n_clusters} clusters")
            with self.console.status(f"Processing {n_clusters} clusters") as status:
                for filename1 in os.listdir(gibbs_matrix_folder):
                    if filename1.endswith(f"of{n_clusters}.mat"):
                        # print(filename1)
                        cluster_found.append(filename1)
                        for filename2 in os.listdir(human_reference_folder):
                            if (
                                hla_list is None
                                or self.formate_HLA_DB(filename2) in hla_list
                            ):
                                correlation = self._compute_and_log_correlation(
                                    gibbs_matrix_folder,
                                    human_reference_folder,
                                    filename1,
                                    filename2,
                                )
                                status.update(
                                    status=f"[bold blue] Compute correlation between {filename1} and {filename2} with correlation {correlation:.4f}",
                                    spinner="squish",
                                    spinner_style="yellow",
                                )
                    else:
                        self.console.log(
                            f"Skipping {filename1} as it does not from {n_clusters} clusters"
                        )

        else:
            self.console.log(
                f"Given n_clusters param {n_clusters} is invalid. Proceeding with 'all' clusters"
            )
            with self.console.status("Processing all clusters") as status:
                for filename1 in os.listdir(gibbs_matrix_folder):
                    for filename2 in os.listdir(human_reference_folder):
                        if (
                            hla_list is None
                            or self.formate_HLA_DB(filename2) in hla_list
                        ):
                            correlation = self._compute_and_log_correlation(
                                gibbs_matrix_folder,
                                human_reference_folder,
                                filename1,
                                filename2,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {filename1} and {filename2} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )

        end_time = time.time()
        elapsed_time = end_time - start_time
        self.console.log(
            f"Cluster Search Preprocess completed in {elapsed_time:.2f} seconds."
        )

    def compute_correlations_v2(
        self,
        db: pd.DataFrame,
        gibbs_results: str,
        n_clusters: str,
        output_path: str,
        hla_list: str = None,
        threshold: float = 0.70,
        treshold_img: bool = False
    ) -> None:
        """
        Compute correlations between test and reference Gibbs matrices.

        :param gibbs_folder: Path to test matrices
        :param human_reference_folder: Path to reference matrices
        :param hla_list: List of HLA types to process (optional, default is all available types)
        """
        #Update to self
        self.threshold = threshold
        self.treshold_img = treshold_img
        self.gibbs_results = gibbs_results
        self.kld_df = read_KLD_file(
            os.path.join(self.gibbs_results, 'images', 'gibbs.KLDvsClusters.tab')
            )
        
        gibbs_result_matrix = os.path.join(gibbs_results, "matrices")
        should_process = False
        if os.path.exists(gibbs_result_matrix) and any(".mat" in file for file in os.listdir(gibbs_result_matrix)):
            self.console.log(
                f"Found Gibbs matrices in {gibbs_result_matrix}")
        else:
            raise FileNotFoundError(
                f"No Gibbs matrices found in {gibbs_result_matrix}")
        # print(db.head())
        self._outfolder = self._make_dir(
            output_path, self.generate_unique_random_ids(6)[0]
        )  # self.generate_unique_random_ids(6)[0]

        if hla_list:
            assert isinstance(
                hla_list, list), "HLA types must be provided as a list [Check main]."
            self.console.log(f"Processing specific HLA types: {hla_list}")
        else:
            hla_list = None
            # self.console.log("Processing all available HLA types.")

        start_time = time.time()
        cluster_found = []
        if n_clusters == "all":
            self.console.log("Processing all clusters", style="blue")
            with self.console.status("Processing all clusters") as status:
                for gibbs_f in os.listdir(gibbs_result_matrix):
                    for mat_path in db['matrices_path']:
                        # print(filename1, filename2,"####"*100)
                        # Extract HLA ID from matrix path

                        if (
                            hla_list is not None
                            and self.formate_HLA_DB(str(mat_path).split('/')[0]) in hla_list
                        ):

                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
                                mat_path,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )
                            # print(correlation)
                        if (
                            hla_list is None

                        ):
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
                                mat_path,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )

                        else:
                            self.console.log(
                                f"Skipping {self.formate_HLA_DB(str(mat_path).split('/')[0])} as it is not in the provided HLA list"
                            )

        elif n_clusters == "best_KL":
            self.console.log("Processing for best KL divergence clusters")
            with self.console.status(
                "Processing best KL divergence clusters"
            ) as status:
                for gibbs_f in os.listdir(gibbs_result_matrix):
                    for mat_path in db['matrices_path']:
                        if (
                            hla_list is not None
                            and self.formate_HLA_DB(str(mat_path).split('/')[0]) in hla_list
                        ):
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
                                mat_path,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )
                        if (
                            hla_list is None

                        ):
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
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

        elif n_clusters.isdigit() and 0 < int(n_clusters) <= 6:

            self.console.log(f"Processing for {n_clusters} clusters")
            with self.console.status(f"Processing {n_clusters} clusters") as status:
                for gibbs_f in os.listdir(gibbs_result_matrix):
                    if gibbs_f.endswith(f"of{n_clusters}.mat"):
                        cluster_found.append(gibbs_f)
                        for mat_path in db['matrices_path']:
                            if (
                                hla_list is not None
                                and self.formate_HLA_DB(str(mat_path).split('/')[-1]) in hla_list
                            ):
                                correlation = self._compute_and_log_correlation_V2(
                                    os.path.join(gibbs_result_matrix, gibbs_f),
                                    mat_path,
                                )
                                status.update(
                                    status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                    spinner="squish",
                                    spinner_style="yellow",
                                )
                            if (
                                hla_list is None

                            ):
                                correlation = self._compute_and_log_correlation_V2(
                                    os.path.join(gibbs_result_matrix, gibbs_f),
                                    mat_path,
                                )
                                status.update(
                                    status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[-1]} with correlation {correlation:.4f}",
                                    spinner="squish",
                                    spinner_style="yellow",
                                )

                            else:
                                self.console.log(
                                    f"Skipping ref databse {self.formate_HLA_DB(str(mat_path).split('/')[-1])} as it is not in the provided HLA list {hla_list}"
                                )
        else:
            self.console.log(
                f"Given n_clusters param {n_clusters} is invalid. Proceeding with 'all' clusters"
            )
            with self.console.status("Processing all clusters") as status:
                for gibbs_f in os.listdir(gibbs_result_matrix):
                    cluster_found.append(gibbs_f)
                    for mat_path in db['matrices_path']:
                        # print(filename1, filename2,"####"*100)
                        if (
                            hla_list is not None
                            or self.formate_HLA_DB(str(mat_path).split('/')[0]) in hla_list
                        ):
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
                                mat_path,
                            )
                            status.update(
                                status=f"[bold blue] Compute correlation between {gibbs_f} and {str(mat_path).split('/')[0]} with correlation {correlation:.4f}",
                                spinner="squish",
                                spinner_style="yellow",
                            )
                        if (
                            hla_list is None

                        ):
                            correlation = self._compute_and_log_correlation_V2(
                                os.path.join(gibbs_result_matrix, gibbs_f),
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
        end_time = time.time()
        elapsed_time = end_time - start_time
        self.console.log(
            f"Cluster Search Preprocess completed in {elapsed_time:.2f} seconds."
        )

        if len(cluster_found) == 0:
            self.console.log(
                f"No cluster files found for {n_clusters} clusters exiting({cluster_found})"
            )
            return

    def _compute_and_log_correlation(
        self,
        gibbs_matrix_folder: str,
        human_reference_folder: str,
        filename1: str,
        filename2: str,
    ) -> None:
        """
        Compute and log the correlation between two Gibbs matrices (test and reference).

        :param gibbs_matrix_folder: Path to test matrices folder
        :param human_reference_folder: Path to reference matrices folder
        :param filename1: File name of the test Gibbs matrix
        :param filename2: File name of the reference Gibbs matrix
        """
        try:
            # Construct full paths to the files
            file1 = os.path.join(gibbs_matrix_folder, filename1)
            file2 = os.path.join(human_reference_folder, filename2)
            # breakpoint()
            # Format input data
            m1 = self.format_input_gibbs(file1)
            m2 = self.format_input_gibbs(file2)

            # Align amino acid order
            m1, m2 = self.amino_acid_order_identical(m1, m2)

            # Calculate correlation
            correlation = m1.corrwith(m2, axis=1).mean()

            # Log the correlation
            # logging.info(f"Correlation between {filename1} and {filename2}: {correlation:.4f}")
            # Store the result in the correlation dictionary
            self.correlation_dict[(filename1, filename2)] = correlation
            return correlation

        except Exception as e:
            # logging.error(f"Failed to compute correlation between {filename1} and {filename2}: {str(e)}")
            self.console.print(
                f"Failed to compute correlation between {filename1} and {filename2}: {str(e)}"
            )
            return float('nan')

    def _compute_and_log_correlation_V2(
        self,
        gibbs_f: str,
        ref_mat: str,
    ) -> None:
        """
        Compute and log the correlation between two Gibbs matrices (test and reference).

        :param gibbs_matrix_folder: Path to test matrices folder
        :param human_reference_folder: Path to reference matrices folder
        :param filename1: File name of the test Gibbs matrix
        :param filename2: File name of the reference Gibbs matrix
        """
        try:
            if not os.path.exists(os.path.join(self._outfolder, 'corr-data')):
                os.makedirs(os.path.join(self._outfolder, 'corr-data'))
            # breakpoint()
            # Format input data
            m1 = self.format_input_gibbs(gibbs_f)
            m2 = self.format_input_gibbs(os.path.join(self.data_dir, ref_mat))

            # Align amino acid order
            m1, m2 = self.amino_acid_order_identical(m1, m2)

            # print(m1, m2)

            # Calculate correlation
            correlation = m1.corrwith(m2, axis=1).mean()
            # print(correlation)
            # breakpoint()

            # Log the correlation
            # logging.info(f"Correlation between {filename1} and {filename2}: {correlation:.4f}")
            # Store the result in the correlation dictionary
            # self.correlation_dict[(str(gibbs_f).split(
            #     '/')[-1], str(ref_mat).split('/')[-1])] = correlation

            self.correlation_dict[(gibbs_f), (ref_mat)] = correlation

            return correlation

        except Exception as e:
            # logging.error(f"Failed to compute correlation between {filename1} and {filename2}: {str(e)}")
            self.console.print(
                f"Failed to compute correlation between {gibbs_f} and {ref_mat}: {str(e)}"
            )
            return float('nan')

    def _find_highest_correlation(self) -> tuple[str, str, float]:
        """
        Find the highest correlation in the correlation dictionary.

        :return: Tuple containing the highest correlation pair and the correlation value
        """
        max_correlation = max(self.correlation_dict.values())
        max_correlation_pair = max(
            self.correlation_dict, key=self.correlation_dict.get)

        return max_correlation_pair[0], max_correlation_pair[1], max_correlation

    def plot_heatmap(self, output_path) -> None:
        """
        Plot a heatmap for the computed correlation dictionary.
        """
        # Check if correlation_dict is empty
        if not self.correlation_dict:

            raise ValueError(
                "correlation_dict is empty. Cannot generate heatmap.[Tip: try running without --n_clusters flag or check your gibbs output files]")

        rows = sorted(set(key[0] for key in self.correlation_dict.keys()))
        cols = sorted(set(key[1] for key in self.correlation_dict.keys()))

        # Create an empty DataFrame with the given rows and columns
        matrix = pd.DataFrame(index=rows, columns=cols, dtype=float)
        # matrix.to_csv(os.path.join(self._outfolder, 'corr-data', 'corr_matrix.csv'),index=False)
        # Fill the matrix with the correlation values
        for (row, col), value in self.correlation_dict.items():
            matrix.loc[row, col] = value

        # Handle NaN values by filling them with 0 or a suitable value
        matrix.fillna(0, inplace=True)

        # Check if the matrix is empty after filling NaN values
        if matrix.empty or matrix.isnull().all().all():
            raise ValueError(
                "Matrix is empty or filled with NaN values. Cannot generate heatmap.(This for calcualting Hihest correaltion)[Tip: try running without --n_clusters flag]")

        custom_cmap = LinearSegmentedColormap.from_list(
            "CustomColours",
            [
                "white", "white", "white", "white", "white", "white", "whitesmoke",
                "lightgrey", "grey", "deeppink",
            ],
        )

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

        ax.vlines(HLA_A_count, y_min, y_max, color="k")
        ax.vlines(HLA_A_count + HLA_B_count, y_min, y_max, color="k")

        unique_names = matrix.index.unique()
        base_names = unique_names.str.replace(r"_gibbs\..*", "", regex=True)
        occurrence_counts = {base: unique_names.str.startswith(
            base).sum() for base in base_names}

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

        # Finalize and save the heatmap
        plt.tight_layout()
        plt.savefig(f"{self._outfolder}/heatmap.png")
        plt.close()

    ##### From here on, Logo comapre module #####

    def find_highest_correlation_for_each_row(self, correlation_dict) -> dict:
        """
        Find the highest correlation for each row in the matrix, i.e., the column with the highest correlation
        for each input sample (row) to the reference HLA type (column).

        :return: A dictionary where each key is the row (cluster/sample) and the value is a tuple
                containing the column (HLA reference) with the highest correlation and the correlation value.
        """
        rows = sorted(set(key[0] for key in correlation_dict.keys()))
        cols = sorted(set(key[1] for key in correlation_dict.keys()))

        # Create a DataFrame to store correlation values
        matrix = pd.DataFrame(index=rows, columns=cols, dtype=float)

        # Populate the correlation matrix
        for (row, col), value in correlation_dict.items():
            matrix.loc[row, col] = value

        highest_corr_per_row = {}

        # Iterate over each row and find the column with the highest correlation
        for row in matrix.index:
            highest_col = matrix.loc[
                row
            ].idxmax()  # Find column with highest correlation for the row
            highest_corr = matrix.loc[
                row, highest_col
            ]  # Get the highest correlation value for the row

            highest_corr_per_row[row] = (highest_col, highest_corr)
            # logging.info(f"Highest correlation for {row} is with {highest_col} with a value of {highest_corr}")

        return highest_corr_per_row

    def add_title(self, img, title, position=(0, 0), font_size=50):
        """
        Add a title to the image at the specified position with an optional font size.
        """
        draw = ImageDraw.Draw(img)

        try:
            font = ImageFont.truetype("arial.ttf", font_size)
        except IOError:
            font = ImageFont.load_default()

        draw.text(position, title, fill="black", font=font)

        return img

    def _find_gibbs_image_path(self, matrix_name, image_folder):
        """
        Find all image paths in the specified folder.
        """
        # print("#####"*100)
        # print(matrix_name)
        for filename in os.listdir(image_folder):
            if filename.endswith(".png") and matrix_name in filename:
                return os.path.join(image_folder, filename)
            elif filename.endswith(".jpg") and matrix_name in filename:
                return os.path.join(image_folder, filename)
            # else:
            #     logging.warning(f"No image found for matrix {matrix_name}")
        return None

    def formate_HLA_DB(self, HLA: str) -> str:
        """
        Format the HLA database file name to extract meaningful parts.

        :param HLA: HLA file name
        :return: Formatted HLA type
        """
        return (
            HLA.replace("HLA_", "").replace(
                "*", "").replace("txt", "").replace(".", "")
        )

    def formate_HLA_user_in(self, HLA: str) -> str:
        """
        Format the HLA database file name to extract meaningful parts.

        :param HLA: HLA file name
        :return: Formatted HLA type
        """
        assert isinstance(HLA, str), "HLA type must be a string."

        if isinstance(HLA, str):
            HLA = HLA.split(",")
            HLA = [
                x.replace("HLA_", "")
                .replace("*", "")
                .replace("txt", "")
                .replace(".", "")
                for x in HLA
            ]
            return HLA
        else:
            return (
                HLA.replace("HLA_", "")
                .replace("*", "")
                .replace("txt", "")
                .replace(".", "")
            )

    def _naturally_presented_log(self, HLA_name, DB_image_folder):
        """
        Find all image paths in the specified folder.
        """
        for filename in os.listdir(DB_image_folder):
            if filename.endswith(".png") and HLA_name in filename:
                return os.path.join(DB_image_folder, filename)
        logging.warning(f"No image found for HLA name {HLA_name}")
        return None

    def _naturally_presented_log_V2(self, HLA, path):
        """
        Find all image paths in the specified folder.
        """
        if os.path.exists(os.path.join(path, f'{HLA}.png')):
            return os.path.join(path, f'{HLA}.png')
        return None

    def create_image_grid(
        self,
        correlation_dict,
        image_folder,
        DB_image_folder,
        output_path,
        columns=3,
        HLA_list=[],
    ):
        """
        Create an image grid where each row corresponds to comparing two images side by side.

        :param correlation_dict: Dictionary of correlations
        :param image_folder: Path to the folder containing the images
        :param DB_image_folder: Path to the folder containing the HLA images
        :param output_path: Path to save the final image grid
        :param columns: Number of columns in the image grid (default is 2 for side-by-side comparison)
        :param HLA_list: List of HLA names provided by the user
        """
        # Get the highest correlation for each row (cluster)
        highest_corr_per_row = self.find_highest_correlation_for_each_row(
            correlation_dict
        )
        # print(highest_corr_per_row)
        display_search_results(highest_corr_per_row, self.threshold)

        sorted_HLA_list = sorted(
            HLA_list,
            key=lambda x: (
                x[0],
                int(re.sub(r"\D", "", x)),
            ),  # Sort by the first letter and then by numeric part
        )
        # Prepare images
        images = []
        titles = []
        highest_corr_per_row_sorted = sorted(
            highest_corr_per_row.items(),
            key=lambda x: x[1][0]
            .split("_")[1]
            # Extract the part after "HLA_" and before ".txt"
            .replace(".txt", ""),
        )

        # print(highest_corr_per_row_sorted)
        for row, (col, corr) in highest_corr_per_row_sorted:
            # Generate title based on row and column
            title = f"{row} -> {col}: {corr:.2f}"
            # title2 = f"naturally presented logo of{self.formate_HLA_DB(col)} -> {row}: {corr:.2f}"
            title2 = f"naturally presented logo of {self.formate_HLA_DB(col)} -> {row}: {corr:.2f}"

            # Find the image paths
            image_path = self._find_gibbs_image_path(
                row.split(".")[1], image_folder)
            nat_img = self._naturally_presented_log(
                self.formate_HLA_DB(col), DB_image_folder
            )
            logging.info(
                f"Cluster matrix {row.split('.')[1]} best allotype match {self.formate_HLA_DB(col)} with correlation {corr:.2f}"
            )

            # Check if the column HLA is in the user provided list, if not find another image or empty
            if self.formate_HLA_DB(col.replace(".", "")) not in HLA_list:
                # Find a user provided HLA alternative if available
                u_hla = [
                    u_hla
                    for u_hla in HLA_list
                    if u_hla.startswith(self.formate_HLA_DB(col)[0])
                ]
                if u_hla:
                    u_hla_nat_img = self._naturally_presented_log(
                        self.formate_HLA_DB(u_hla[0]), DB_image_folder
                    )
                    title3 = (
                        f"{self.formate_HLA_DB(u_hla[0])} -> {row}: {corr:.2f}"
                        if u_hla
                        else "Provided HLA type found"
                    )
                else:
                    u_hla_nat_img = None
                    title3 = "Provided HLA type found"
            else:
                u_hla_nat_img = None
                title3 = "Provided HLA type found"

            # Open the images
            if image_path and nat_img:
                img1 = Image.open(image_path)
                img2 = Image.open(nat_img)
                # title3 = f"{formate_HLA_DB(u_hla[0])} -> {row}: {corr:.2f}" if u_hla else "Provided HLA type found"

                if u_hla_nat_img:
                    img3 = Image.open(u_hla_nat_img)
                else:
                    img3 = Image.new(
                        "RGB", (img1.width, img1.height), color=(255, 255, 255)
                    )  # Empty image

                img1 = self.add_title(
                    img1, title, position=(img1.width // 4, 5), font_size=60
                )  # Increased font size
                img2 = self.add_title(
                    img2, title2, position=(img2.width // 4, 5), font_size=60
                )  # Increased font size
                img3 = self.add_title(
                    img3, title3, position=(img2.width // 4, 5), font_size=60
                )  # Increased font size

                # Append both images
                images.append(img1)
                images.append(img2)
                images.append(img3)
        if not images:
            self.console.log(
                "Error: No images were found to create a grid.", style="bold red")
            return
        # Calculate grid dimensions
        rows = (len(images) + columns - 1) // columns
        width = images[0].width * columns
        height = images[0].height * rows

        grid_image = Image.new("RGB", (width, height))

        # Paste the images into the grid, 2 images per row
        for i, img in enumerate(images):
            row = i // columns
            col = i % columns
            grid_image.paste(img, (col * img.width, row * img.height))

        # Save the final image grid
        grid_image.save(f"{self._outfolder}/compare_allotypes.png")
        # logging.info(f"Output saved in {self._outfolder}")
        self.console.log(f"Output saved in {self._outfolder}")

    # New Image Grid Module
    def generate_image_grid(self, correlation_dict):
        """
        Generate an image grid for the correlation results.
        """
        highest_corr_per_row = self.find_highest_correlation_for_each_row(
            correlation_dict
        )
        # print(highest_corr_per_row)
        display_search_results(highest_corr_per_row, self.threshold)
        # print(highest_corr_per_row)
        # config_dict = {
        #     "image_folder": "data/ref_data/Gibbs_motifs_human/logos",
        #     "DB_image_folder": "data/ref_data/Gibbs_motifs_human/images",
        #     "output_path": "data/output",
        #     "columns": 3,
        #     "HLA_list": [],
        # }

        # imagelayout()
    def _make_correlation_plot(self, gibbs_df, motif_df, mat_motif):
        df1 = pd.DataFrame(gibbs_df)
        df2 = pd.DataFrame(motif_df)

        df1['Position'] = df1.index + 1
        df2['Position'] = df2.index + 1

        df1_melted = df1.melt(id_vars='Position',
                              var_name='Amino Acid', value_name='Cluster')
        df2_melted = df2.melt(id_vars='Position',
                              var_name='Amino Acid', value_name='Reference')

        # Merge the two dataframes
        df_merged = pd.merge(df1_melted, df2_melted, on=[
                             'Position', 'Amino Acid'])

        # Create the base chart
        base = alt.Chart(df_merged, width="container").mark_circle().encode(
            x='Cluster',
            y='Reference',
            color='Amino Acid',
            tooltip=['Amino Acid', 'Cluster', 'Reference', 'Position']
        )
        # Calculate the correlation coefficient
        corr_coef = df_merged[['Cluster', 'Reference']].corr().iloc[0, 1]
        # corr_coef = df_merged['Cluster'].corrwith(df_merged['Reference']).iloc[0]
        # corr_coef = pd.Series(df_merged['Cluster']).corr(pd.Series(df_merged['Reference']))

        # corrwith

        # Create the regression line
        regression_line = base.transform_regression('Cluster', 'Reference').mark_line(opacity=0.50, shape='mark').transform_fold(
            ["reg-line"],
            as_=["Regression", "y"]
        ).encode(alt.Color("Regression:N"))

        # Add the correlation coefficient as text
        corr_text = alt.Chart(pd.DataFrame({
            'Cluster': [df_merged['Cluster'].min()],
            'Reference': [df_merged['Reference'].max()],
            'text': [f'Correlation: {corr_coef:.2f}']
        })).mark_text(align='left', baseline='top', dx=8, dy=-5).encode(
            x='Cluster:Q',
            y='Reference:Q',
            text='text:N'
        )
        # Combine the charts
        chart = base + regression_line + corr_text
        if not os.path.exists(os.path.join(self._outfolder, 'corr-data')):
            os.makedirs(os.path.join(self._outfolder, 'corr-data'))
        chart.save(
            f"{os.path.join(self._outfolder, 'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.json")
        chart.save(
            f"{os.path.join(self._outfolder, 'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.png")

    def insert_script_hla_section(self, script_data_path, div_id):
        script_template = Template('''
        fetch('{{ script_data_path }}')
            .then(response => response.json())
            .then(data => {
            console.log(data);
            // Process the data as needed
            var opt = {"renderer": "canvas", "actions": false};
            vegaEmbed("#{{ div_id }}", data, opt);
            })
            .catch(error => console.error('Error fetching the JSON data:', error));
        ''')
        return script_template.render(script_data_path=script_data_path, div_id=div_id)

    def insert_script_png_json(self, script_data_path, img_fallback_path, div_id):
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

        return script_template.render(script_data_path=script_data_path, img_fallback_path=img_fallback_path, div_id=div_id)

    def render_hla_section(self, hla_name, corr, best_cluster_img, naturally_presented_img,kld_clust_group_kld):
        
        if str(self.species).lower() == "human":
            hla_name = f"HLA-{hla_name}"
        else:
            hla_name
        # if kld == None:
        #     kld = "None"
        template = Template('''
        <div class="row" style="border: 2px solid #007bff;">
        <div class="row">
            <h5 style="text-align: center;">Best Matched Allotype is {{ hla_name }} and PCC = {{corr}} , KLD= {{ kld }}</h5>
        </div>
        <div class="row">
            <div class="col">
            <div class="card">
                <h5 style="text-align: center;">Gibbs Cluster</h5>
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
                <h5 style="text-align: center;">Reference Motif logo of {{ hla_name }}</h5>
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
        return template.render(hla_name=hla_name, corr=corr, best_cluster_img=best_cluster_img, naturally_presented_img=naturally_presented_img,kld=kld_clust_group_kld)

    def make_datatable(self, correlation_dict):
        df = pd.DataFrame(correlation_dict.items(),
                          columns=['HLA', 'Correlation'])
        df['Cluster'] = df['HLA'].apply(
            lambda x: x[0].split('/')[-1].split('.')[1])
        df['HLA'] = df['HLA'].apply(lambda x: x[1])
        df['HLA'] = df['HLA'].apply(
            lambda x: x.split('/')[-1].replace('.txt', ''))
        df['Correlation'] = df['Correlation'].apply(lambda x: round(x, 2))
        #Add KLD from kld_clust_group_kld
        if self.kld_df is not None:
            try:
                df['KLD'] = df['Cluster'].apply(
                    lambda x: self.kld_df.loc[
                        self.kld_df['cluster'] == int(str(x).split('of')[-1]),
                        f'group{str(x).split("of")[0]}'
                    ].values[0] if f'group{str(x).split("of")[0]}' in self.kld_df.columns else None
                )
            except KeyError as e:
                # Handle the case where the key is not found
                self.console.log(
                    f"KeyError: The key KLD group was not found in the Corr DataFrame."
                )
                df['KLD'] = 'NA'
        else:
            df['KLD'] = 'NA'
            

        df = df.sort_values(by='Correlation', ascending=False)
        df = df.reset_index(drop=True)
        df = df[['Cluster', 'HLA', 'Correlation', 'KLD']]
        
        # self.console.log(
        #     f"Correlation data:\n{df}", style="bold green")
        fig, ax = plot_kld_pcc_distance(self.kld_df,df)
        fig.savefig(os.path.join(self._outfolder, 'corr-data', 'kld_pcc_distance.png'))
        self.console.log(
            f"Saved KLD and PCC distance plot at {os.path.join(self._outfolder, 'corr-data', 'kld_pcc_distance.png')}", style="bold green")
        return df

    def process_correlation_data(self, df=None):
        """
        Process the correlation matrix data from a CSV file and generate the necessary JSON for visualization
        """
        if df is not None:
            # self.console.log(f"Processing correlation data from {df}")
            # Extract unique clusters and HLA types
            unique_clusters = sorted(df['Cluster'].unique())
            unique_hlas = sorted(df['HLA'].unique())

            self.console.log(
                f"Found {len(unique_clusters)} unique clusters and {len(unique_hlas)} unique HLA types", style="bold green")

            # Create a pivot table for the heatmap
            pivot_df = df.pivot_table(
                index='Cluster',
                columns='HLA',
                values='Correlation',
                fill_value=0  # Fill missing values with 0
            )

            # Create heatmap data for D3.js visualization
            heatmap_data = []
            for cluster in unique_clusters:
                for hla in unique_hlas:
                    # Get correlation value if it exists
                    try:
                        correlation = df.loc[(df['Cluster'] == cluster) & (
                            df['HLA'] == hla), 'Correlation'].values[0]
                    except IndexError:
                        correlation = 0  # Default if no correlation exists

                    heatmap_data.append({
                        "cluster": cluster,
                        "hla": hla,
                        "correlation": correlation
                    })

            # Create a matrix format for alternative visualization
            matrix_data = []
            for cluster in unique_clusters:
                cluster_data = {"cluster": cluster, "values": []}
                for hla in unique_hlas:
                    try:
                        correlation = df.loc[(df['Cluster'] == cluster) & (
                            df['HLA'] == hla), 'Correlation'].values[0]
                    except IndexError:
                        correlation = 0  # Default if no correlation exists
                    cluster_data["values"].append(correlation)
                matrix_data.append(cluster_data)

            # Create visualization data object
            visualization_data = {
                "heatmapData": heatmap_data,
                "matrixData": matrix_data,
                "clusters": unique_clusters,
                "hlas": unique_hlas
            }

            # Save the data as JSON
            output_file = Path("hla_correlation_data.json")
            with open(os.path.join(self._outfolder, output_file), 'w') as f:
                json.dump(visualization_data, f, indent=2)

            self.console.log(f"Correlation data saved as {output_file}")

            # Create a static heatmap image as a reference
            # self.create_static_heatmap(pivot_df)

            return visualization_data

        else:
            self.console.log("No correlation data found to process")
            return None

    def save_correlation_data(self):
        """Save correlation data to CSV and JSON files"""
        try:
            # Save dataframe to CSV
            self.corr_df.to_csv(os.path.join(
                self._outfolder, 'corr-data', 'corr_matrix.csv'), index=False)
            self.console.log(f"Saved PCC data at corr-data/corr_matrix.csv",style="bold green")
            # Save to JSON format for D3.js
            try:
                json_data = self.write_to_json(self.corr_df)
                with open(os.path.join(self._outfolder, 'corr-data', 'pcc_data.json'), 'w+') as f:
                    json.dump(json_data, f, ensure_ascii=False, indent=4)
                self.console.log(f"Saved Network and Heatmap data at corr-data/pcc_data.json" ,style="bold green")
                return True

            except Exception as e:
                self.console.log(
                    f"Failed to save the correlation matrix JSON for network plot: {str(e)}",style="bold red")
                return False
                 
        except Exception as e:
            self.console.log(
                f"Failed to save the correlation matrix: {str(e)}",style="bold red")
            return False

    def write_to_json(self, df):
        """Convert dataframe to JSON format for D3.js visualization"""
        data = []

        for idx, row in df.iterrows():
            # Convert each row to the required format
            data.append({
                "cluster_id": row['Cluster'],
                "hla_type": row['HLA'],
                "correlation": float(row['Correlation'])
            })

        return data

    def make_heatmap(self, correlation_dict):
        df = self.make_datatable(correlation_dict)
        # print(df)

        try:
            self.corr_df = df
            self.d3js_json = self.save_correlation_data()
            return True
            # self.d3js_json = self.process_correlation_data(df)
        except Exception as e:
            self.console.log(
                f"Failed to save the save_correlation_data() -- fun : {str(e)}")
        return False

    def make_datatable_html(self, correlation_dict, df=None):
        if df is None:
            df = self.make_datatable(correlation_dict)
            table_start = """
                <div class="container py-4">
        <div class="card">
            <div class="card-header bg-secondary text-white">
                <h4 class="card-title mb-0">Correlation For all Reference HLA/MHC</h4>
            </div>
            <div class="card-body">
            
            <table id="correlation_table"  class="table table-bordered">
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
                if df['Correlation'][i] >= self.threshold:
                    tr += f"""
                    <tr class="table-success">
                        <td>{i+1}</td>
                        <td>{df['Cluster'][i]}</td>
                        <td>{df['HLA'][i]}</td>
                        <td>{df['Correlation'][i]}</td>
                        <td>{round(df['KLD'][i], 4) if df['KLD'][i] != 'NA' else 'NA'}</td>
                    </tr>
                    """
                if df['Correlation'][i] >= 0.5 and df['Correlation'][i] < self.threshold:
                    tr += f"""
                    <tr class="table-warning">
                        <td>{i+1}</td>
                        <td>{df['Cluster'][i]}</td>
                        <td>{df['HLA'][i]}</td>
                        <td>{df['Correlation'][i]}</td>
                        <td>{round(df['KLD'][i], 4) if df['KLD'][i] != 'NA' else 'NA'}</td>
                    </tr>
                    """
                if df['Correlation'][i] < 0.5:
                    tr += f"""
                    <tr class="table-danger">
                        <td>{i+1}</td>
                        <td>{df['Cluster'][i]}</td>
                        <td>{df['HLA'][i]}</td>
                        <td>{df['Correlation'][i]}</td>
                        <td>{round(df['KLD'][i], 4) if df['KLD'][i] != 'NA' else 'NA'}</td>
                    </tr>
                    """

            table_end = """
         
               </tbody>
            </table>
            
               </div>
             </div>
              </div>
            """
        return table_start + tr + table_end

        # return df.to_html(classes='table table-striped', index=False, table_id='correlation_table')

    #!!! Best HLA View
    
    def generate_best_matches_tab(self, highest_corr_per_row, threshold=0.7):
        """
        Extracts the best HLA-A, HLA-B, and HLA-C matches from highest_corr_per_row and generates the HTML tab view.
        
        Args:
            highest_corr_per_row: Dictionary mapping cluster ID to tuple of (hla_id, correlation)
            threshold: Minimum correlation value to consider (default: 0.7)
            
        Returns:
            HTML string for the Best Matched Results tab content
        """
        # Initialize dictionaries to store best matches for each HLA type
        best_hla_matches = {
            'A': {'best': None, 'second': None},
            'B': {'best': None, 'second': None},
            'C': {'best': None, 'second': None}
        }
        
        # Parse highest_corr_per_row to extract HLA type information
        for cluster_id, (hla_id, correlation) in highest_corr_per_row.items():
            # Extract HLA type (A, B, or C) from the HLA ID
            # Check different formats: HLA_A0101, HLA-A0101, or A0101
            hla_name = hla_id
            if isinstance(hla_id, str):
                if 'HLA_' in hla_id:
                    hla_name = hla_id.split('HLA_')[-1]
                elif 'HLA-' in hla_id:
                    hla_name = hla_id.split('HLA-')[-1]
                
                # Extract the HLA type letter (A, B, C)
                hla_type = None
                if hla_name.startswith('A'):
                    hla_type = 'A'
                elif hla_name.startswith('B'):
                    hla_type = 'B'
                elif hla_name.startswith('C'):
                    hla_type = 'C'
                
                # Skip if correlation is below threshold or HLA type not recognized
                if correlation < threshold or hla_type is None:
                    continue
                    
                # Format HLA name for display
                formatted_hla = f"HLA-{hla_name}"
                
                # Check if this is better than current best match
                if best_hla_matches[hla_type]['best'] is None or correlation > best_hla_matches[hla_type]['best'][1]:
                    # Move current best to second best
                    best_hla_matches[hla_type]['second'] = best_hla_matches[hla_type]['best']
                    # Set new best
                    best_hla_matches[hla_type]['best'] = (formatted_hla, correlation, cluster_id)
                # Check if this is better than current second best
                elif best_hla_matches[hla_type]['second'] is None or correlation > best_hla_matches[hla_type]['second'][1]:
                    best_hla_matches[hla_type]['second'] = (formatted_hla, correlation, cluster_id)
        
        # Generate HTML for tab content
        html = '<div class="tab-pane fade show active" id="beastresults" role="tabpanel" aria-labelledby="beastresults-tab">'
        
        # Function to lookup KLD value (to be implemented based on your data structure)
        def get_kld_value(cluster_id, hla_id):
            # Replace with actual lookup in your KLD dataframe or data structure
            # For now returning a placeholder value
            return 1.5
        
        # Function to get image paths (to be implemented based on your data structure)
        def get_image_paths(cluster_id, hla_id):
            # Replace with actual path construction logic
            gibbs_img = f"cluster-img/gibbs_logos_{str(cluster_id).replace(".mat","")}-001.png"
            ref_img = f"allotypes-img/{str(hla_id).replace(".txt","").replace("-","_")}.png"
            return gibbs_img, ref_img
        
        # Generate row for each HLA type
        for hla_type, color in [('A', 'primary'), ('B', 'success'), ('C', 'info')]:
            if best_hla_matches[hla_type]['best'] is None:
                continue  # Skip if no matches found for this HLA type
            
            html += f'<!-- HLA-{hla_type} Row -->\n<div class="row">'
            
            # First column (best match)
            if best_hla_matches[hla_type]['best']:
                hla_id, correlation, cluster_id = best_hla_matches[hla_type]['best']
                kld = get_kld_value(cluster_id, hla_id)
                gibbs_img, ref_img = get_image_paths(cluster_id, hla_id)
                
                html += f'''
                <div class="col-md-6">
                    <div class="card mt-4 shadow-sm">
                        <div class="card-header bg-{color} text-white">
                            <h5 class="card-title mb-0 text-center">{hla_id} (PCC: {correlation:.2f}, KLD: {kld:.2f})</h5>
                        </div>
                        <div class="card-body">
                            <div class="row mb-3">
                                <h5 class="mb-2">Gibbs Cluster</h5>
                                <div class="text-center">
                                    <img src="{gibbs_img}" class="img-fluid border" alt="Gibbs Cluster" style="max-height: 300px;">
                                </div>
                            </div>
                            <div class="row">
                                <h5 class="mb-2">Reference HLA</h5>
                                <div class="text-center">
                                    <img src="{ref_img}" class="img-fluid border" alt="Reference HLA" style="max-height: 300px;">
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                '''
            
            # Second column (second best match)
            if best_hla_matches[hla_type]['second']:
                hla_id, correlation, cluster_id = best_hla_matches[hla_type]['second']
                kld = get_kld_value(cluster_id, hla_id)
                gibbs_img, ref_img = get_image_paths(cluster_id, hla_id)
                
                html += f'''
                <div class="col-md-6">
                    <div class="card mt-4 shadow-sm">
                        <div class="card-header bg-{color} text-white">
                            <h5 class="card-title mb-0 text-center">{hla_id} (PCC: {correlation:.2f}, KLD: {kld:.2f})</h5>
                        </div>
                        <div class="card-body">
                            <div class="row mb-3">
                                <h5 class="mb-2">Gibbs Cluster</h5>
                                <div class="text-center">
                                    <img src="{gibbs_img}" class="img-fluid border" alt="Gibbs Cluster" style="max-height: 300px;">
                                </div>
                            </div>
                            <div class="row">
                                <h5 class="mb-2">Reference HLA</h5>
                                <div class="text-center">
                                    <img src="{ref_img}" class="img-fluid border" alt="Reference HLA" style="max-height: 300px;">
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                '''
            
            html += '</div>\n'
        
        html += '</div>'
        
        return html

    def create_tab_views(self, highest_corr_per_row, all_results_html):
        """
        Creates the complete tabbed interface with Best Matched Results and All Results tabs.
        
        Args:
            highest_corr_per_row: Dictionary mapping cluster ID to tuple of (hla_id, correlation)
            all_results_html: HTML content for the All Results by Cluster tab
            
        Returns:
            String containing the complete HTML for the tabbed interface
        """
        # Generate the best matches tab content
        best_matches_html = self.generate_best_matches_tab(highest_corr_per_row)
        
        # Create the tabs structure
        tabs_html = f'''
        
        <!-- Tabs navigation -->
        <ul class="nav nav-tabs" id="compareTabs" role="tablist">
            <li class="nav-item" role="presentation">
                <button class="nav-link active" id="beastresults-tab" data-bs-toggle="tab"
                    data-bs-target="#beastresults" type="button" role="tab" aria-controls="combined"
                    aria-selected="true">Best Matched Results</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link active" id="gibbsresults-tab" data-bs-toggle="tab"
                    data-bs-target="#gibbsresults" type="button" role="tab" aria-controls="combined"
                    aria-selected="false">Gibbs KLD</button>
            </li>
            <li class="nav-item" role="presentation">
                <button class="nav-link" id="allresults-tab" data-bs-toggle="tab"
                    data-bs-target="#allresults" type="button" role="tab" aria-controls="heatmap"
                    aria-selected="false">All Results by Cluster</button>
            </li>
        </ul>
        
        <div class="tab-content" id="compareTabsContent">
            <!-- Best Matched Results Tab -->
            {best_matches_html}
            
            <!-- Gibbs KLD Tab -->
            <div class="tab-pane fade" id="gibbsresults" role="tabpanel" aria-labelledby="gibbsresults-tab">
            <div class="col-md-12">
                    <div class="card mt-4 shadow-sm">
                        <div class="card-body">
                            <div class="row">
                                <div class="text-center">
                                    <img src="corr-data/kld_pcc_distance.png" class="img-fluid border" alt="Gibbs Cluster" style="max-height: 400px;">
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            <!-- All Results by Cluster Tab -->
            <div class="tab-pane fade" id="allresults" role="tabpanel" aria-labelledby="allresults-tab">
                {all_results_html}
            </div>
        </div>
        '''
        
        return tabs_html


    #### NEW !!!!! Carousel Control based clusters #########


    def create_cluster_hierarchy(self, highest_corr_per_row, gibbs_out):
        """
        Create a hierarchical dictionary of clusters organized by their groups.
        """
        result = {}
        
        for row_path, (col_path, correlation) in highest_corr_per_row.items():
            try:
                filename = os.path.basename(row_path)
                
                if 'gibbs.' not in filename:
                    continue
                
                cluster_id = filename.split('gibbs.')[1].split('.')[0]
                group_num, cluster_num = map(int, cluster_id.split('of'))
                
                hla = os.path.basename(col_path).replace('.txt', '').split('_')[1] if self.species == 'human' else os.path.basename(col_path).replace('.txt', '')
                
                if cluster_num not in result:
                    result[cluster_num] = {}
                
                result[cluster_num][group_num] = {
                    'id': cluster_id,
                    'matrix_path': row_path,
                    'motif': col_path,
                    'hla': hla,
                    'correlation': float(correlation)
                }
                
                gibbs_img = self._find_gibbs_image_path(cluster_id, os.path.join(gibbs_out, 'logos'))
                if gibbs_img:
                    dest_path = os.path.join(self._outfolder, 'cluster-img', os.path.basename(gibbs_img))
                    shutil.copy(gibbs_img, dest_path)
                    result[cluster_num][group_num]['gibbs_img'] = dest_path
                
                if hasattr(self, 'db_path') and hasattr(self, 'db'):
                    db_matches = self.db[self.db['formatted_allotypes'] == hla]
                    if not db_matches.empty:
                        nat_path = db_matches['motif_path'].values[0]
                        if nat_path:
                            src_path = os.path.join(self.db_path, nat_path)
                            dest_path = os.path.join(self._outfolder, 'allotypes-img', os.path.basename(nat_path))
                            shutil.copy(src_path, dest_path)
                            result[cluster_num][group_num]['nat_img'] = dest_path
            
            except Exception as e:
                if hasattr(self, 'console'):
                    self.console.print(f"Error processing cluster {row_path}: {str(e)}")
        
        return result

    def render_cluster_carousels(self, highest_corr_per_row, gibbs_out):
        """
        Renders carousels for clusters, grouping them by their group number.
        """
        cluster_hierarchy = self.create_cluster_hierarchy(highest_corr_per_row, gibbs_out)
        # cluster_html = """
        #         <div class="container py-4">
        #         <div class="card">
        #             <div class="card-header bg-secondary text-white">
        #                 <h4 class="card-title mb-0">Gibbs Cluster and Reference Motif Comparisons</h4>
        #             </div>
        #             <div class="card-body">
        # """
        
        cluster_html = """ 
                    """
        
    
        
        
        for cluster_num in sorted(cluster_hierarchy.keys()):
            if cluster_hierarchy[cluster_num]:
                cluster_html += self._create_carousel_for_cluster(
                    f"carousel-cluster-{cluster_num}", cluster_num, cluster_hierarchy[cluster_num])

        if self.treshold_img:
            cluster_html += f"""
                 </div>
                 <h2> Note: Filtered cluster groups below threshold {self.threshold}</h2>
                 </div>
                        </div>
                """
        else:
            cluster_html += """
                        </div>
                        </div>
                        </div>
                """
        
        return cluster_html


    def _create_carousel_for_cluster(self, carousel_id, cluster_num, group_data):
        """
        Creates a Bootstrap carousel for a specific cluster group.

        Args:
            carousel_id: Unique ID for the carousel
            cluster_num: The cluster number
            group_data: Dictionary of groups in this cluster

        Returns:
            HTML string for the carousel
        """
        
        # self.treshold_img = True
        # Filter groups based on threshold if treshold_img is enabled
        if self.treshold_img:
            group_data = {k: v for k, v in group_data.items() if v.get('correlation', 0) >= self.threshold}
            
        # Get ordered list of group numbers
        group_nums = sorted(group_data.keys())

        # Skip if no groups to display
        if not group_nums:
            return ""

        total_slides = len(group_nums)

        # Start building the carousel HTML
        carousel_html = f"""
        <div class="row mt-4 mb-4">
            <div class="col-12">
            <h2 class="text-center">RESULTS for {cluster_num} CLUSTERS ({'1 group' if str(cluster_num) == "1" else f'1 to {cluster_num} groups'})</h2>
            <div id="{carousel_id}" class="carousel slide" data-ride="carousel" data-interval="false">
        """

        # Add indicators
        carousel_html += """
                    <ol class="carousel-indicators">
        """

        for i, group_num in enumerate(group_nums):
            active_class = "active" if i == 0 else ""
            carousel_html += f'<li data-bs-target="#{carousel_id}" data-bs-slide-to="{i}" class="{active_class}"></li>\n'

        carousel_html += """
                    </ol>
                    <div class="carousel-inner">
        """
        # Add KLD score if available
        
        
        if self.kld_df is not None:
            kld_clust_df = self.kld_df[self.kld_df['cluster'] == cluster_num]
            if not kld_clust_df.empty:
                kld = kld_clust_df['total'].values[0]
                self.console.log(f"-------KLD Results for {cluster_num} clusters-------", style="bold green")
                self.console.log(f"Total KLD score for cluster {cluster_num}: {kld}", style="bold yellow")
            else:
                kld_clust_df = None
        
        # Add carousel items
        for i, group_num in enumerate(group_nums):
            active_class = "active" if i == 0 else ""
            cluster_data = group_data[group_num]

            # Get HLA info
            hla_name = cluster_data.get('hla', 'Unknown HLA')
            correlation = cluster_data.get('correlation', 0)
            correlation_formatted = round(correlation, 2) if isinstance(
                correlation, (int, float)) else correlation

            # Get image paths
            gibbs_img = cluster_data.get('gibbs_img', None)
            nat_img = cluster_data.get('nat_img', None)

            # Remove base path if needed
            if gibbs_img and hasattr(self, '_outfolder'):
                gibbs_img = str(gibbs_img).replace(f"{self._outfolder}/", '')

            if nat_img and hasattr(self, '_outfolder'):
                nat_img = str(nat_img).replace(f"{self._outfolder}/", '')

            #find KLD for group 
            kld_clust_group_kld = kld_clust_df[f'group{group_num}'].values[0] if kld_clust_df is not None and f'group{group_num}' in kld_clust_df.columns else None
            self.console.log(f"Group {group_num} KLD score for cluster {cluster_num}: {kld_clust_group_kld}", style="bold yellow")
            # Generate the HLA section for this cluster
            kld_clust_group_kld = round(kld_clust_group_kld, 2) if isinstance(
                kld_clust_group_kld, (int, float)) else kld_clust_group_kld
            
            hla_section = self.render_hla_section(
                hla_name, correlation_formatted, gibbs_img, nat_img,kld_clust_group_kld)

            # Create carousel item with the HLA section
            carousel_html += f"""
                        <div class="carousel-item {active_class}">
                            <div class="container">
                                {hla_section}
                                <div class="text-center mt-2">Slide {i+1} of {total_slides}</div>
                            </div>
                        </div>
            """

        # Add carousel controls
        carousel_html += f"""
                    </div>
                    <a class="carousel-control-prev" href="#{carousel_id}" role="button" data-slide="prev">
                        <span class="carousel-control-prev-icon" aria-hidden="true"></span>
                        <span class="sr-only">Previous</span>
                    </a>
                    <a class="carousel-control-next" href="#{carousel_id}" role="button" data-slide="next">
                        <span class="carousel-control-next-icon" aria-hidden="true"></span>
                        <span class="sr-only">Next</span>
                    </a>
                </div>
            </div>
        </div>
        """

        return carousel_html


    def render_clustered_results_nav(self, highest_corr_per_row, gibbs_out):
        """
        Renders the complete HTML for all clustered results with carousels and necessary JavaScript.
        """
        
        html_card = """
                <div class="container py-4">
                <div class="card">
                    <div class="card-header bg-secondary text-white">
                        <h4 class="card-title mb-0">Gibbs Cluster and Reference Motif Comparisons</h4>
                    </div>
                    <div class="card-body">
                    
        
        """
        
        # cluster_hierarchy = self.create_cluster_hierarchy(highest_corr_per_row, gibbs_out)

        all_results_html = self.render_cluster_carousels(highest_corr_per_row, gibbs_out)
        
        
        # Add the carousel HTML to the card
        html_card += self.create_tab_views(
            highest_corr_per_row, all_results_html)
        
        html_card += """"
                </div>
                </div>
                </div>
        """
        
        return html_card







    # New functions ends here

    def generate_html_layout(self, correlation_dict, db, gibbs_out, immunolyser=False):
        """
        Generate an image grid for the correlation results.
        """
        highest_corr_per_row = self.find_highest_correlation_for_each_row(
            correlation_dict
        )
        # print(highest_corr_per_row)
        display_search_results(highest_corr_per_row, self.threshold)
        # print(highest_corr_per_row)

        # added to self
        self.db = db

        if not os.path.exists(os.path.join(self._outfolder, 'cluster-img')):
            os.makedirs(os.path.join(self._outfolder, 'cluster-img'))
        if not os.path.exists(os.path.join(self._outfolder, 'allotypes-img')):
            os.makedirs(os.path.join(self._outfolder, 'allotypes-img'))
        output_dict = {
        }
        # D3 heatmap
        heatmap_d3_html = """
    <div class="container py-4">
        <div class="card">
            <div class="card-header bg-secondary text-white">
                <h4 class="card-title mb-0">HLA Correlation Analysis</h4>
            </div>
            <div class="card-body">
                <!-- Threshold control -->
                <div class="mb-4">
                    <label for="threshold" class="form-label">Correlation Threshold: <span id="threshold-value">0.70</span></label>
                    <input type="range" class="form-range" id="threshold" min="0" max="1" step="0.01" value="0.70">
                </div>
                
                <!-- Tabs navigation -->
                <ul class="nav nav-tabs" id="visualizationTabs" role="tablist">
                    <li class="nav-item" role="presentation">
                        <button class="nav-link active" id="combined-tab" data-bs-toggle="tab" data-bs-target="#combined" type="button" role="tab" aria-controls="combined" aria-selected="true">Combined View</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="heatmap-tab" data-bs-toggle="tab" data-bs-target="#heatmap" type="button" role="tab" aria-controls="heatmap" aria-selected="false">Heatmap</button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link" id="network-tab" data-bs-toggle="tab" data-bs-target="#network" type="button" role="tab" aria-controls="network" aria-selected="false">Network</button>
                    </li>
                </ul>
                
                <!-- Tab content -->
                <div class="tab-content" id="visualizationTabsContent">
                    <div class="tab-pane fade show active" id="combined" role="tabpanel" aria-labelledby="combined-tab">
                        <div id="heatmap-container-combined" class="visualization-container" style="min-height: 300px; position: relative;"></div>
                        <div id="network-container-combined" class="visualization-container" style="min-height: 400px; position: relative;">
                            <button class="btn btn-sm btn-secondary center-network-btn" onclick="centerNetwork('network-container-combined')">Center View</button>
                        </div>
                    </div>
                    <div class="tab-pane fade" id="heatmap" role="tabpanel" aria-labelledby="heatmap-tab">
                        <div id="heatmap-container-single" class="visualization-container" style="min-height: 500px; position: relative;"></div>
                    </div>
                    <div class="tab-pane fade" id="network" role="tabpanel" aria-labelledby="network-tab">
                        <div id="network-container-single" class="visualization-container" style="min-height: 600px; position: relative;">
                            <button class="btn btn-sm btn-secondary center-network-btn" onclick="centerNetwork('network-container-single')">Center View</button>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    """

        html_create = html_content + body_start

        carousel_js = """
        $(document).ready(function() {
            // Initialize all carousels
            $('.carousel').carousel({
                interval: false  // Prevent auto-sliding
            });
            
            // Fix for carousel controls
            $('.carousel-control-prev, .carousel-control-next').click(function(e) {
                e.preventDefault();
                var targetId = $(this).attr('href');
                var direction = $(this).hasClass('carousel-control-prev') ? 'prev' : 'next';
                $(targetId).carousel(direction);
            });
            
            // Fix for indicators
            $('.carousel-indicators li').click(function() {
                var targetId = $(this).data('target');
                var slideIndex = $(this).data('slide-to');
                $(targetId).carousel(parseInt(slideIndex));
            });
        });
        """
        body_end_1 = Template("""
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
  <script src ="https://cdn.datatables.net/2.2.2/js/dataTables.bootstrap5.js"></script>

<footer class="text-body-secondary py-5">
    <div class="container">
      <p class="float-end mb-1">
        <a href="#">Back to top</a>
      </p>
      <p class="mb-1"> Clust-Search v{{version}}</p>
      <p class="mb-1"> Li Lab / Purcell Lab, Monash University&copy; 2025,please refere git repo for more information</p>
      <p class="mb-0">github<a href="https://github.com/Sanpme66/HLA-PepClust">Visit git hub the homepage</a> or read our <a
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


""").render(version = __version__)

        br_tag = """
        <br>
"""

        body_end_2 = """
</script>


    
</body>
</html>
"""
        
        heatmap_d3_js = """
        
        
        // Load PCC data from JSON file
let data = [];

// Fetch and process JSON data
fetch('corr-data/pcc_data.json')
    .then(response => response.json())
    .then(jsonData => {
        // Process the data into the format expected by the visualization
        data = jsonData.map(item => ({
            Cluster: item.cluster_id,
            HLA: item.hla_type,
            Correlation: item.correlation
        }));
        
        // Now that data is loaded, initialize visualizations
        
        const simulations = {};

        // Color mapping for clusters
        const clusterColors = {
            "1of1": "#FFFF33",
            "1of2": "#FF7F00",
            "2of2": "#FF7F00",
            "1of3": "#984EA3",
            "2of3": "#984EA3",
            "3of3": "#984EA3",
            "1of4": "#4DAF4A",
            "2of4": "#4DAF4A",
            "3of4": "#4DAF4A",
            "4of4": "#4DAF4A",
            "1of5": "#377EB8",
            "2of5": "#377EB8",
            "3of5": "#377EB8",
            "4of5": "#377EB8",
            "5of5": "#377EB8",
            "1of6": "#E41A1C",
            "2of6": "#E41A1C",
            "3of6": "#E41A1C",
            "4of6": "#E41A1C",
            "5of6": "#E41A1C",
            "6of6": "#E41A1C",
        };
        
        // Initialize tooltip
        const tooltip = d3.select("body")
            .append("div")
            .attr("class", "tooltip-custom")
            .style("opacity", 0);
        
        // Current threshold value
        let threshold = 0.70;
        
        // Initialize visualizations
        createHeatmap("heatmap-container-combined", threshold);
        createNetwork("network-container-combined", threshold);
        createHeatmap("heatmap-container-single", threshold);
        createNetwork("network-container-single", threshold);
        
        // Set up threshold slider
        const thresholdSlider = document.getElementById("threshold");
        const thresholdValueDisplay = document.getElementById("threshold-value");
        
        thresholdSlider.addEventListener("input", function() {
            threshold = parseFloat(this.value);
            thresholdValueDisplay.textContent = threshold.toFixed(2);
            
            // Update visualizations
            createHeatmap("heatmap-container-combined", threshold);
            createNetwork("network-container-combined", threshold);
            
            // Only update the visible tab's visualization if not on combined view
            const activeTab = document.querySelector(".tab-pane.active");
            if (activeTab.id === "heatmap") {
                createHeatmap("heatmap-container-single", threshold);
            } else if (activeTab.id === "network") {
                createNetwork("network-container-single", threshold);
            }
        });
        
        // Update visualizations when switching tabs
        const tabs = document.querySelectorAll('button[data-bs-toggle="tab"]');
        tabs.forEach(tab => {
            tab.addEventListener('shown.bs.tab', function(event) {
                const targetId = event.target.getAttribute("aria-controls");
                if (targetId === "heatmap") {
                    createHeatmap("heatmap-container-single", threshold);
                } else if (targetId === "network") {
                    createNetwork("network-container-single", threshold);
                }
            });
        });
        
        // Create heatmap visualization
        function createHeatmap(containerId, threshold) {
            const container = document.getElementById(containerId);
            
            // Clear previous content
            container.innerHTML = '';
            
            // Filter data
            const filteredData = data.filter(d => d.Correlation >= threshold);
            
            if (filteredData.length === 0) {
                const noDataMsg = document.createElement("div");
                noDataMsg.className = "text-center p-4";
                noDataMsg.textContent = "No data to display with current threshold";
                container.appendChild(noDataMsg);
                return;
            }
            
            // Get unique clusters and HLAs
            const clusters = [...new Set(filteredData.map(d => d.Cluster))].sort();
            const hlaTypes = [...new Set(filteredData.map(d => d.HLA))].sort();
            
            // Set up dimensions
            const containerWidth = container.clientWidth || 800;
            
            // Calculate responsive cell size based on available width
            const maxCellSize = 40;
            const minCellSize = 25;
            const availableWidth = containerWidth - 150; // Account for margins and labels
            const calculatedCellSize = Math.max(minCellSize, Math.min(maxCellSize, availableWidth / hlaTypes.length));
            const cellSize = calculatedCellSize;
            
            // Increased bottom margin to accommodate labels, and adjusted right margin
            const margin = { top: 50, right: 20, bottom: 180, left: 100 };
            
            // Calculate dimensions with adjusted margins
            const width = cellSize * hlaTypes.length + margin.left + margin.right;
            const height = cellSize * clusters.length + margin.top + margin.bottom;
            
            // Create a div with overflow visible for the SVG
            const svgContainer = document.createElement('div');
            svgContainer.className = 'svg-container';
            container.appendChild(svgContainer);
            
            // Create SVG
            const svg = d3.select(svgContainer)
                .append("svg")
                .attr("width", width)
                .attr("height", height)
                .attr("viewBox", `0 0 ${width} ${height}`)
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);
            
            // Create scales
            const x = d3.scaleBand()
                .domain(hlaTypes)
                .range([0, cellSize * hlaTypes.length])
                .padding(0.05);
            
            const y = d3.scaleBand()
                .domain(clusters)
                .range([0, cellSize * clusters.length])
                .padding(0.05);
            
            // Color scale
            const colorScale = d3.scaleSequential(d3.interpolateYlGnBu)
                .domain([0, 1]);
            
            // Add Y axis labels
            svg.append("g")
                .selectAll("text")
                .data(clusters)
                .enter()
                .append("text")
                .attr("x", -10)
                .attr("y", d => y(d) + y.bandwidth() / 2)
                .attr("text-anchor", "end")
                .attr("dominant-baseline", "middle")
                .text(d => d)
                .style("font-size", "12px")
                .style("font-weight", "500");
            
            // Create labels group for X axis
            const labelGroup = svg.append("g").attr("class", "axis-labels");
            
            // Compact one-line labels with vertically staggered positions
            hlaTypes.forEach((d, i) => {
                const xPos = x(d) + x.bandwidth() / 2;
                const hlaType = d.replace("HLA_", "");
                const isEven = i % 2 === 0;
                
                // Single line compact label
                labelGroup.append("text")
                    .attr("x", xPos)
                    .attr("y", isEven ? -10 : -30)
                    .attr("text-anchor", "middle")
                    .attr("dominant-baseline", "middle")
                    .text(hlaType)
                    .style("font-size", "10px")
                    .style("font-weight", "500");
                
                // Add connecting line
                labelGroup.append("line")
                    .attr("x1", xPos)
                    .attr("y1", isEven ? -5 : -25)
                    .attr("x2", xPos)
                    .attr("y2", -5)
                    .style("stroke", "#ddd")
                    .style("stroke-width", 0.5)
                    .style("stroke-dasharray", "2,2");
            });
            
            // Create cells
            svg.selectAll("rect")
                .data(filteredData)
                .enter()
                .append("rect")
                .attr("x", d => x(d.HLA))
                .attr("y", d => y(d.Cluster))
                .attr("width", x.bandwidth())
                .attr("height", y.bandwidth())
                .style("fill", d => colorScale(d.Correlation))
                .style("stroke", "#fff")
                .style("stroke-width", "1px")
                .on("mouseover", function(event, d) {
                    // Highlight cell
                    d3.select(this)
                        .style("stroke", "#333")
                        .style("stroke-width", "2px");
                    
                    // Show tooltip
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", 1);
                        
                    tooltip.html(`
                        <strong>Cluster:</strong> ${d.Cluster}<br>
                        <strong>HLA:</strong> ${d.HLA}<br>
                        <strong>Correlation:</strong> ${d.Correlation.toFixed(2)}
                    `)
                    .style("left", (event.pageX + 10) + "px")
                    .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    // Remove highlight
                    d3.select(this)
                        .style("stroke", "#fff")
                        .style("stroke-width", "1px");
                    
                    // Hide tooltip
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                });
            
            // Add cell text
            svg.selectAll("text.cell-text")
                .data(filteredData)
                .enter()
                .append("text")
                .attr("class", "cell-text")
                .attr("x", d => x(d.HLA) + x.bandwidth() / 2)
                .attr("y", d => y(d.Cluster) + y.bandwidth() / 2)
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "middle")
                .text(d => d.Correlation.toFixed(2))
                .style("fill", d => d.Correlation > 0.5 ? "white" : "black")
                .style("font-size", "11px")
                .style("font-weight", "bold");
        }
        
        // Create network visualization
        function createNetwork(containerId, threshold) {
            const container = document.getElementById(containerId);
            
            // Clear previous content
            container.innerHTML = '';
            
            // Add center button
            const centerBtn = document.createElement("button");
            centerBtn.className = "btn btn-sm btn-secondary center-network-btn";
            centerBtn.textContent = "Center View";
            centerBtn.onclick = function() { centerNetwork(containerId); };
            container.appendChild(centerBtn);
            
            // Filter data
            const filteredData = data.filter(d => d.Correlation >= threshold);
            
            if (filteredData.length === 0) {
                const noDataMsg = document.createElement("div");
                noDataMsg.className = "text-center p-4";
                noDataMsg.textContent = "No data to display with current threshold";
                container.appendChild(noDataMsg);
                return;
            }
            
            // Process data for network
            const clusters = [...new Set(filteredData.map(d => d.Cluster))];
            const hlaTypes = [...new Set(filteredData.map(d => d.HLA))];
            
            // Create nodes
            const nodes = [
                ...clusters.map(c => ({ id: c, group: "cluster" })),
                ...hlaTypes.map(h => ({ id: h, group: "hla" }))
            ];
            
            // Create links
            const links = filteredData.map(d => ({
                source: d.Cluster,
                target: d.HLA,
                value: d.Correlation
            }));
            
            // Set dimensions
            const width = container.clientWidth || 800;
            const height = container.clientHeight || 400;
            
            // Create SVG
            const svg = d3.select(container)
                .append("svg")
                .attr("width", width)
                .attr("height", height);
            
            // Add zoom behavior
            const zoom = d3.zoom()
                .scaleExtent([0.5, 5])
                .on("zoom", (event) => {
                    g.attr("transform", event.transform);
                });
            
            svg.call(zoom);
            
            // Create a container for the graph that can be transformed
            const g = svg.append("g");
            
            // Create simulation
            const simulation = d3.forceSimulation(nodes)
                .force("link", d3.forceLink(links).id(d => d.id).distance(d => 100 - d.value * 40))
                .force("charge", d3.forceManyBody().strength(-100))
                .force("center", d3.forceCenter(width / 2, height / 2))
                .force("collision", d3.forceCollide().radius(d => d.group === "cluster" ? 25 : 10))
                .force("x", d3.forceX(width / 2).strength(0.05))
                .force("y", d3.forceY(height / 2).strength(0.05));
            
            // Store simulation reference
            simulations[containerId] = simulation;
            
            // Create links
            const link = g.append("g")
                .selectAll("line")
                .data(links)
                .enter()
                .append("line")
                .attr("stroke-width", d => Math.max(1, d.value * 3))
                .attr("stroke", d => clusterColors[d.source.id] || "#999")
                .attr("opacity", 0.6);
            
            // Create nodes
            const node = g.append("g")
                .selectAll("circle")
                .data(nodes)
                .enter()
                .append("circle")
                .attr("r", d => d.group === "cluster" ? 15 : 6)
                .attr("fill", d => {
                    if (d.group === "cluster") {
                        return clusterColors[d.id] || "#999";
                    }
                    return "#ccc";
                })
                .attr("stroke", "#fff")
                .attr("stroke-width", 1.5)
                .on("mouseover", function(event, d) {
                    // Highlight node
                    d3.select(this)
                        .attr("stroke", "#333")
                        .attr("stroke-width", 2);
                    
                    // Show tooltip
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", 1);
                    
                    if (d.group === "cluster") {
                        // Show cluster info
                        const clusterData = filteredData.filter(item => item.Cluster === d.id);
                        const avgCorrelation = clusterData.reduce((sum, item) => sum + item.Correlation, 0) / clusterData.length;
                        
                        tooltip.html(`
                            <strong>${d.id}</strong><br>
                            Average Correlation: ${avgCorrelation.toFixed(3)}<br>
                            HLA Connections: ${clusterData.length}
                        `);
                    } else {
                        // Show HLA connections
                        const connections = filteredData.filter(item => item.HLA === d.id)
                            .sort((a, b) => b.Correlation - a.Correlation);
                        
                        let html = `<strong>${d.id}</strong><br><table style="border-collapse: collapse;">
                            <tr><th style="padding: 3px; text-align: left;">Cluster</th>
                            <th style="padding: 3px; text-align: right;">Correlation</th></tr>`;
                        
                        connections.forEach(conn => {
                            html += `<tr>
                                <td style="padding: 3px;">${conn.Cluster}</td>
                                <td style="padding: 3px; text-align: right;">${conn.Correlation.toFixed(2)}</td>
                            </tr>`;
                        });
                        
                        html += `</table>`;
                        tooltip.html(html);
                    }
                    
                    tooltip.style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                        
                    // Highlight connections
                    if (d.group === "cluster") {
                        link.attr("opacity", l => l.source.id === d.id ? 1 : 0.1)
                            .attr("stroke-width", l => l.source.id === d.id ? Math.max(2, l.value * 4) : 1);
                    } else {
                        link.attr("opacity", l => l.target.id === d.id ? 1 : 0.1)
                            .attr("stroke-width", l => l.target.id === d.id ? Math.max(2, l.value * 4) : 1);
                    }
                })
                .on("mouseout", function() {
                    // Remove highlight
                    d3.select(this)
                        .attr("stroke", "#fff")
                        .attr("stroke-width", 1.5);
                    
                    // Hide tooltip
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                    
                    // Reset link highlighting
                    link.attr("opacity", 0.6)
                        .attr("stroke-width", d => Math.max(1, d.value * 3));
                })
                .call(d3.drag()
                    .on("start", dragstarted)
                    .on("drag", dragged)
                    .on("end", dragended)
                );
            
            // Add labels to cluster nodes
            g.append("g")
                .selectAll("text")
                .data(nodes.filter(d => d.group === "cluster"))
                .enter()
                .append("text")
                .attr("text-anchor", "middle")
                .attr("dominant-baseline", "middle")
                .text(d => d.id)
                .style("font-size", "12px")
                .style("fill", "white")
                .style("pointer-events", "none");
            
            // Add labels to HLA nodes with high correlation
            g.append("g")
                .selectAll("text.hla-label")
                .data(nodes.filter(d => {
                    if (d.group !== "hla") return false;
                    
                    // Find max correlation for this HLA
                    const maxCorr = Math.max(...filteredData
                        .filter(item => item.HLA === d.id)
                        .map(item => item.Correlation));
                    
                    return maxCorr > 0.7;
                }))
                .enter()
                .append("text")
                .attr("class", "hla-label")
                .attr("dx", 8)
                .attr("dy", 3)
                .text(d => d.id.replace("HLA_", ""))
                .style("font-size", "8px")
                .style("fill", "#333")
                .style("pointer-events", "none");
            
            // Update positions on simulation tick
            simulation.on("tick", () => {
                // Constrain nodes to container boundaries
                nodes.forEach(d => {
                    // Add padding based on node size
                    const padding = d.group === "cluster" ? 20 : 10;
                    
                    d.x = Math.max(padding, Math.min(width - padding, d.x));
                    d.y = Math.max(padding, Math.min(height - padding, d.y));
                });
                
                link
                    .attr("x1", d => d.source.x)
                    .attr("y1", d => d.source.y)
                    .attr("x2", d => d.target.x)
                    .attr("y2", d => d.target.y);
                
                node
                    .attr("cx", d => d.x)
                    .attr("cy", d => d.y);
                
                g.selectAll("text")
                    .filter(function() {
                        return !this.classList.contains("hla-label");
                    })
                    .data(nodes.filter(d => d.group === "cluster"))
                    .attr("x", d => d.x)
                    .attr("y", d => d.y);
                
                g.selectAll("text.hla-label")
                    .data(nodes.filter(d => d.group === "hla" && filteredData
                        .filter(item => item.HLA === d.id)
                        .some(item => item.Correlation > 0.7)))
                    .attr("x", d => d.x)
                    .attr("y", d => d.y);
            });
            
            // Position nodes in a more balanced layout initially
            setTimeout(() => {
                // Start with cluster nodes in a circle
                const clusterNodes = nodes.filter(d => d.group === "cluster");
                const radius = Math.min(width, height) / 3;
                const angleStep = (2 * Math.PI) / clusterNodes.length;
                
                clusterNodes.forEach((d, i) => {
                    const angle = i * angleStep;
                    d.x = width / 2 + radius * Math.cos(angle);
                    d.y = height / 2 + radius * Math.sin(angle);
                    d.fx = d.x;
                    d.fy = d.y;
                });
                
                // Run simulation for a short time
                simulation.alpha(1).restart();
                
                // Release fixed positions after brief arrangement
                setTimeout(() => {
                    clusterNodes.forEach(d => {
                        d.fx = null;
                        d.fy = null;
                    });
                    simulation.alpha(0.3).restart();
                }, 1000);
            }, 100);
            
            // Drag functions
            function dragstarted(event, d) {
                if (!event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            }
            
            function dragged(event, d) {
                d.fx = event.x;
                d.fy = event.y;
            }
            
            function dragended(event, d) {
                if (!event.active) simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            }
        }
        
        // Center network view function
        function centerNetwork(containerId) {
            const container = document.getElementById(containerId);
            const svg = d3.select(container).select("svg");
            const simulation = simulations[containerId];
            
            if (svg.empty() || !simulation) return;
            
            // Reset zoom
            svg.transition().duration(750).call(
                d3.zoom().transform,
                d3.zoomIdentity.translate(0, 0).scale(1)
            );
            
            // Reset node positions
            const width = container.clientWidth;
            const height = container.clientHeight;
            
            simulation.force("center", d3.forceCenter(width / 2, height / 2))
                .alpha(1)
                .restart();
        }
    })
    .catch(error => {
        console.error('Error loading correlation data:', error);
        document.querySelectorAll('.visualization-container').forEach(container => {
            container.innerHTML = `
                <div class="alert alert-danger" role="alert">
                    Error loading data: json_data_path not found or invalid format.\n
                    Tip: Please use VS code live server it loads pcc_data.json
                </div>
            `;
        });
    });
    
"""
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

        # heatmap_json = self.make_heatmap(correlation_dict)
        self.make_heatmap(correlation_dict)


        # script_js = ""
        df_corr_html = self.make_datatable_html(correlation_dict)
        # print(df_corr_html)
        # print(df_co)
        immunolyser_out = ""
        immunolyser_out_js = ""
        # cluster = {}
        # print("##"*100)
        # print(highest_corr_per_row)

        # new_html_carousel = self.render_cluster_carousels(
        #     highest_corr_per_row, gibbs_out)
        new_html_carousel = self.render_clustered_results_nav(
            highest_corr_per_row, gibbs_out)
        body_end_1 += carousel_js
        # print(new_html_carousel)

        # print("$$$"*100)

        # for row, (col, corr) in highest_corr_per_row.items():
        #     # cluster = str(row).split("/")[-1].split('.')[1]
        #     if corr >= self.threshold:
        #         output_dict[str(row).split("/")[-1].split('.')[1]] = {
        #             'cluster': str(row).split("/")[-1],
        #             'HLA': str(col).split('/')[-1].replace('.txt',''),
        #             'correlation': corr,
        #             'gibbs_img': None,
        #             'nat_img': None,
        #             'corr_plot': None,
        #             'corr_json': None
        #         }
        #         gibbs_img = self._find_gibbs_image_path(str(row).split("/")[-1].split('.')[1], os.path.join(gibbs_out, 'logos'))

        #         if gibbs_img:
        #             # shutil.copy(gibbs_img, os.path.join(self._outfolder, 'cluster-img', f"{str(gibbs_img).split('/')[-1]}"))
        #             # Use os.path functions for platform-independent path handling
        #             gibbs_img_basename = os.path.basename(gibbs_img)
        #             dest_path = os.path.join(os.path.join(self._outfolder, 'cluster-img'), gibbs_img_basename)
        #             shutil.copy(gibbs_img, dest_path)

        #             # output_dict[str(row).split("/")[-1].split('.')[1]]['gibbs_img'] = os.path.join(self._outfolder, 'cluster-img', f"{str(gibbs_img).split('/')[-1]}")
        #             # Use platform-independent path operations for the output dictionary
        #             row_basename = os.path.basename(str(row))
        #             row_name = os.path.splitext(row_basename)[0].split('.')[1]
        #             output_dict[row_name]['gibbs_img'] = dest_path

        #         if self.species == 'human':
        #             hla = str(col).split('/')[-1].replace('.txt','').split('_')[1]
        #         else:
        #             hla = str(col).split('/')[-1].replace('.txt','')
        #         nat_path = db[db['formatted_allotypes'] == hla]['motif_path'].values[0]

        #         if nat_path:
        #             # shutil.copy(os.path.join(self.db_path,nat_path), os.path.join(self._outfolder, 'allotypes-img', f"{str(nat_path).split('/')[-1]}"))
        #             src_path = os.path.join(self.db_path, nat_path)
        #             nat_path_basename = os.path.basename(nat_path)
        #             dest_path = os.path.join(os.path.join(self._outfolder, 'allotypes-img'), nat_path_basename)
        #             # Copy the file
        #             shutil.copy(src_path, dest_path)
        #             # output_dict[str(row).split("/")[-1].split('.')[1]]['nat_img'] = os.path.join(self._outfolder, 'allotypes-img', f"{str(nat_path).split('/')[-1]}")
        #             row_basename = os.path.basename(str(row))
        #             row_name = os.path.splitext(row_basename)[0].split('.')[1]
        #             output_dict[row_name]['nat_img'] = dest_path

        #         try:
        #             gibbs_mt = self.format_input_gibbs(row)
        #             nat_mat = self.format_input_gibbs(os.path.join(self.db_path,col))

        #             # Align amino acid order
        #             gibbs_mt, nat_mat = self.amino_acid_order_identical(gibbs_mt, nat_mat)
        #             mat_motif = str(row).split("/")[-1].split('.')[1]+"_"+str(col).split('/')[-1].replace('.txt','')
        #             # Remove Corelation plot
        #             # self._make_correlation_plot(gibbs_mt, nat_mat,mat_motif)
        #             # output_dict[str(row).split("/")[-1].split('.')[1]]['corr_plot'] = f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.png"
        #             # output_dict[str(row).split("/")[-1].split('.')[1]]['corr_json'] = f"{os.path.join(self._outfolder,'corr-data')}/amino_acids_comparison_with_correlation_{mat_motif}.json"

        #         except Exception as e:
        #             self.console.print(
        #                 f"Failed to compute correlation plot between {row} and {col}: {str(e)}"
        #             )
        #             # return float('nan')
        #         rows_list = self.render_hla_section(hla,round(corr,2), str(output_dict[str(row).split("/")[-1].split('.')[1]]['gibbs_img']).replace(f"{self._outfolder}/",''), str(output_dict[str(row).split("/")[-1].split('.')[1]]['nat_img']).replace(f"{self._outfolder}/",''))
        #         html_create += rows_list
        #         # plot_js = self.insert_script_hla_section(str(output_dict[str(row).split("/")[-1].split('.')[1]]['corr_json']).replace(f"{self._outfolder}/",''), f"correlation_chart_{str(row).split('/')[-1].split('.')[1]}")
        #         plot_js = self.insert_script_png_json(str(output_dict[str(row).split("/")[-1].split('.')[1]]['corr_json']).replace(f"{self._outfolder}/",''),str(output_dict[str(row).split("/")[-1].split('.')[1]]['corr_plot']).replace(f"{self._outfolder}/",''),f"correlation_chart_{str(row).split('/')[-1].split('.')[1]}")
        #         body_end_1  += plot_js

        #         immunolyser_out += rows_list
        #         immunolyser_out_js += plot_js

        # if heatmap_json:
        #     body_end_1 += heatmap_js
        #     immunolyser_out_js += heatmap_js
        # else:
        #     heatmap_div = ""
        #     heatmap_js =""
        if self.d3js_json:
            body_end_1 += heatmap_d3_js
            immunolyser_out_js += heatmap_d3_js
        else:
            heatmap_d3_html = """
        <div class="container py-4">
            <div class="card">
                <div class="card-header bg-secondary text-white">
                    <h4 class="card-title mb-0">HLA Correlation Analysis</h4>
                </div>
                <div class="card-body">
                <h2> Faild to load Heatmap and Network plot</h2>
                </div>
            </div>
        </div>
            """
            heatmap_d3_js = ""
        # addeing new new_html_carousel
        html_create += new_html_carousel
        # ends here
        html_create += br_tag + \
            br_tag + heatmap_d3_html + \
            br_tag + df_corr_html + br_tag + body_end_1 + body_end_2

        immunolyser_out += br_tag + br_tag + df_corr_html
        with open(os.path.join(self._outfolder, "clust-search-result.html"), "w") as file:
            file.write(html_create)

        if immunolyser:
            with open(os.path.join(self._outfolder, "immunolyser-out.html"), "w") as file:
                file.write(immunolyser_out)
            with open(os.path.join(self._outfolder, "immunolyser-out.js"), "w") as file:
                file.write(immunolyser_out_js)
        self.console.log(
            f"HTML layout saved in {os.path.join(self._outfolder, 'clust-search-result.html')}", style="bold green")

        # image_path = self._find_gibbs_image_path(row.split(".")[1], image_folder)
        # nat_img = self._naturally_presented_log(self.formate_HLA_DB(col), DB_image_folder)
        #     if gibbs_img and nat_path:
        #         html_content += f"""
        #         <div>
        #             <img src="{output_dict[str(row).split("/")[-1].split('.')[1]]['gibbs_img']}" alt="{row}">
        #             <p>{row} -> {col}: {corr:.2f}</p>
        #         </div>
        #         <div>
        #             <img src="{output_dict[str(row).split("/")[-1].split('.')[1]]['nat_img']}" alt="{col}">
        #             <p>naturally presented logo of {self.formate_HLA_DB(col)} -> {row}: {corr:.2f}</p>
        #         </div>
        #         """

        # html_content += """
        #     </div>
        # </body>
        # </html>
        # """

        # with open(os.path.join(self._outfolder, "clust-search-result.html"), "w") as file:
        #     file.write(html_content)
        # self.console.log(f"HTML layout saved in {os.path.join(self._outfolder, 'clust-search-result.html')}")


def generate_heatmap_html(self, correlation_dict):
    df = self.make_datatable(correlation_dict)
    # print(df)


def run_cluster_search(args):
    # if args.output_folder is None:
    #     os.makedirs("output", exist_ok=True)
    #     output_folder = "output"
    cluster_search = ClusterSearch()
    cluster_search.console.rule(
        "[bold red]Stage 1/2: Data processing for correlation matrices."
    )
    # print(args.species)
    if str(args.species).lower() in ["mouse", "human"]:
        CONSOLE.log(
            f"Species provided: [bold yellow]{args.species}", style="bold green")
        CONSOLE.log(
            f"Loading reference databse for [bold yellow]{args.species}", style="bold green")
        db = cluster_search._db_loader(args.reference_folder, args.species)
        CONSOLE.log(f"Reference database loaded successfully.",
                    style="bold green")
    else:
        raise ValueError(
            "Invalid species provided. Please provide a valid species. refer only. [Human or Mouse]")
    if args.threshold:

        CONSOLE.log(
            f"Threshold provided: [bold yellow]{args.threshold}", style="bold green")
        if args.threshold > 1 or args.threshold < 0:
            raise ValueError(
                "Invalid threshold provided. Please provide a valid threshold between 0 and 1")
    else:
        args.threshold = 0.70
        CONSOLE.log(
            f"No threshold provided. Using default threshold of 0.70", style="bold yellow")

    if args.hla_types:
        CONSOLE.log(
            f"HLA/MHC allotypes types provided: [bold yellow]{args.hla_types}", style="bold green")
        u_hla_list = []
        with CONSOLE.status("Checking HLA/MHC in databse") as status:
            status.update(
                status=f"[bold blue] Loading HLA/MHC in databse",
                spinner="squish",
                spinner_style="yellow",
            )
            for u_hla in str(args.hla_types).split(','):
                if u_hla in db['formatted_allotypes'].values:
                    status.update(
                        status=f"[bold blue] HLA/MHC {u_hla} found in databse",
                        spinner="squish",
                        spinner_style="yellow",
                    )
                    u_hla_list.append(u_hla)
                    CONSOLE.log(
                        f"HLA/MHC [bold yellow]{u_hla}[yellow] found in databse", style="bold green")
                else:
                    status.update(
                        status=f"[bold blue] HLA/MHC {u_hla} not found in databse",
                        spinner="squish",
                        spinner_style="yellow",
                    )
                    CONSOLE.log(
                        f"HLA/MHC [bold yellow]{u_hla} not found in databse", style="bold red")
        if len(u_hla_list) > 0:
            CONSOLE.log(
                f"Valid HLA/MHC allotypes types provided: [bold yellow]{u_hla_list}", style="bold green")
            args.hla_types = u_hla_list
    else:
        args.hla_types = None
        CONSOLE.log(
            f"No HLA/MHC allotypes types provided. Using all available HLA types from the reference folder.", style="bold yellow")

    CONSOLE.log(f"calculating compute_correlations.", style="bold blue")
    cluster_search.compute_correlations_v2(
        db,
        args.gibbs_folder,
        args.n_clusters,
        args.output,
        args.hla_types,
        args.threshold,
        args.treshold_img
    )

    # cluster_search.generate_image_grid(cluster_search.correlation_dict,db)
    cluster_search.generate_html_layout(
        cluster_search.correlation_dict, db, args.gibbs_folder, args.immunolyser)

    # breakpoint()

    # cluster_search.compute_correlations(
    #     args.gibbs_folder,
    #     args.reference_folder,
    #     args.n_clusters,
    #     args.output,
    #     args.hla_types,
    # )

    cluster_search.console.rule(
        "[bold red]Stage 2/2: Finding best matching Naturally presented HLA ."
    )
    # if args.output_folder is None:
    cluster_search.plot_heatmap(args.output)

    # cluster_search.console.rule("[bold red]Stage 3/4: Cheking the HLA.")

    # cluster_search.check_HLA_DB(args.hla_types, args.reference_folder)
    # # cluster_search.create_image_grid(cluster_search.correlation_dict, os.path.join(args.gibbs_folder, 'logos'), os.path.join(args.reference_folder, 'images'), os.path.join(args.output_folder, 'image_grid_D90.png'), HLA_list=cluster_search.valid_HLA)

    # cluster_search.console.rule(
    #     "[bold red]Stage 4/4: Finding Best Matched HLA for Each cluster."
    # )
    # cluster_search.create_image_grid(
    #     cluster_search.correlation_dict,
    #     os.path.join(args.gibbs_folder, "logos"),
    #     str(args.reference_folder).replace(
    #         "/output_matrices_human", "").replace("output_matrices", ""),
    #     os.path.join(args.output, "compare_motif_corr.png"),
    #     HLA_list=cluster_search.valid_HLA,
    # )
    # logging.info("Process completed successfully.")
    if args.log:
        log_file_path = os.path.join(
            cluster_search._outfolder, "search_cluster.log")
        save_console_log()
        # raise FileNotFoundError(f"Log file not found: {log_file_path}")

# Remove thie after test


# if __name__ == "__main__":
    #     # ClusterSearch().compute_correlations(
    #     #     "data/9mersonly",
    #     #     "data/ref_data/Gibbs_motifs_mouse/output_matrices",
    #     #     "all",
    #     #     "data/outputM",
    #     # )

    #     # ClusterSearch()._compute_and_log_correlation(
    #     #     "data/9mersonly",
    #     #     "data/ref_data/Gibbs_motifs_human/output_matrices_human",
    #     #     "cluster_1of5.mat",
    #     #     "HLA_A_02_01.txt",
    #     # )
    #     # print(sys.argv)

    #     # print(ClusterSearch()._db_loader("data/ref_data/","mouse"))

    # run_cluster_search(
    #     argparse.Namespace(
    #         credits=False,
    #         gibbs_folder="data/9mersonly",
    #         species="mouse",
    #         hla_types="H2_Db,H2_Dd,H2_Dq,H2_Kb,H2_Kd,H2_Kk",
    #         log=False,
    #         n_clusters="5",
    #         output="data/outputMouseTest",
    #         processes=4,
    #         version=False,
    #         immunolyser=False
    #     )
    # )
