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
from rich.traceback import install
import sys

install(show_locals=True)
from cli.logger import *
from cli.imagegrid import imagelayout


# Configure logging
# logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# logging = CONSOLE
# _logConfig(logSave=True)
class ClusterSearch:
    def __init__(self):
        self.correlation_dict = {}
        self.valid_HLA = []
        self.console = CONSOLE

    def generate_unique_random_ids(self, count: int) -> list:
        """
        Generate a list of unique 6-digit random IDs.

        :param count: Number of unique random IDs to generate
        :return: List of unique 6-digit random IDs
        """
        start = 100000
        end = 999999
        if count > (end - start + 1):
            raise ValueError("Count is larger than the range of unique IDs available.")

        return np.random.choice(
            range(start, end + 1), size=count, replace=False
        ).tolist()
    def _db_loader(self, db_path: str,species:str) -> pd.DataFrame:
        """
        Load the database file.

        :param db_path: Path to the database file
        :return: DataFrame containing the database
        """
        if not os.path.exists(os.path.join(db_path, f'{species}.db')):
            raise FileNotFoundError(f"Database file {db_path} does not exist.")

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

        raise FileNotFoundError(f"No cluster file found for {n_clusters} clusters.")

    def _make_dir(self, path: str,rand_ids:int) -> None:
        """
        Make a directory if it does not exist.

        :param path: Path to the directory
        """
        # if not os.path.exists(os.path.join(path, f'clust_result_{rand_ids}')):
        #     os.makedirs(os.path.join(path, f'clust_result_{rand_ids}'))
            
        if not os.path.exists(os.path.join(path, f'clust_result')):
            os.makedirs(os.path.join(path, f'clust_result'))
        #return os.path.join(path, f'clust_result_{rand_ids}')
        
        self.console.log(f"Output directory created {os.path.join(path, f'clust_result')}")

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

        for HLA in self.formate_HLA_user_in(HLA_list):
            if self.formate_HLA_DB(HLA) not in DB_hla_list:
                # logging.error(f"HLA type {HLA} not found in the reference folder.")
                self.console.log(f"HLA type {HLA} not found in the reference folder.")
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
        df.iloc[:, 0] = df.iloc[:, 0].str.replace(r"^\d+\s\w\s", "", regex=True)
        new_df = pd.DataFrame(df.iloc[1:, 0].str.split(expand=True).values, columns=amino_acids)
        new_df.reset_index(drop=True, inplace=True)
        new_df = new_df.apply(pd.to_numeric, errors='coerce')
        return new_df

    @staticmethod
    def amino_acid_order_identical(df1: pd.DataFrame, df2: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        if list(df1.columns) != list(df2.columns):
            logging.warning("The amino acid column order is different. Reordering columns.")
            df2 = df2[df1.columns]
        return df1, df2
    
    @staticmethod
    def formate_HLA_DB(HLA: str) -> str:
        """
        Format the HLA database file name to extract meaningful parts.

        :param HLA: HLA file name
        :return: Formatted HLA type
        """
        return HLA.replace("HLA_", "").replace("*", "").replace("txt", "")
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
                self.console.print(f"HLA type {HLA} not found in the reference folder.")
                # logging.error(f"Available HLA types: {DB_hla_list}")
                # self.console.print(f"Available HLA types: {DB_hla_list}")
                return False
            else:
                self.valid_HLA.append(HLA)

        return True

    def compute_correlations(
        self,
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
        self._outfolder = self._make_dir(
            output_path, self.generate_unique_random_ids(6)[0]
        ) #self.generate_unique_random_ids(6)[0]

        if hla_list:
            hla_list = hla_list[0].split(
                ","
            )  # Split if comma-separated string is passed
            assert isinstance(hla_list, list), "HLA types must be provided as a list."
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

        start_time = time.time()

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
                        print(filename1)
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

    def _find_highest_correlation(self) -> tuple[str, str, float]:
        """
        Find the highest correlation in the correlation dictionary.

        :return: Tuple containing the highest correlation pair and the correlation value
        """
        max_correlation = max(self.correlation_dict.values())
        max_correlation_pair = max(self.correlation_dict, key=self.correlation_dict.get)
        return max_correlation_pair[0], max_correlation_pair[1], max_correlation

    def plot_heatmap(self, output_path) -> None:
        """
        Plot a heatmap for the computed correlation dictionary.
        """
        # Check if correlation_dict is empty
        if not self.correlation_dict:
            raise ValueError("correlation_dict is empty. Cannot generate heatmap.")

        rows = sorted(set(key[0] for key in self.correlation_dict.keys()))
        cols = sorted(set(key[1] for key in self.correlation_dict.keys()))
        
        # Create an empty DataFrame with the given rows and columns
        matrix = pd.DataFrame(index=rows, columns=cols, dtype=float)

        # Fill the matrix with the correlation values
        for (row, col), value in self.correlation_dict.items():
            matrix.loc[row, col] = value

        # Handle NaN values by filling them with 0 or a suitable value
        matrix.fillna(0, inplace=True)

        # Check if the matrix is empty after filling NaN values
        if matrix.empty or matrix.isnull().all().all():
            raise ValueError("Matrix is empty or filled with NaN values. Cannot generate heatmap.")

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
        occurrence_counts = {base: unique_names.str.startswith(base).sum() for base in base_names}

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
            
            if filename.endswith(".jpg") and matrix_name in filename:
                return os.path.join(image_folder, filename)
        logging.warning(f"No image found for matrix {matrix_name}")
        return None

    def formate_HLA_DB(self, HLA: str) -> str:
        """
        Format the HLA database file name to extract meaningful parts.

        :param HLA: HLA file name
        :return: Formatted HLA type
        """
        return (
            HLA.replace("HLA_", "").replace("*", "").replace("txt", "").replace(".", "")
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
        display_search_results(highest_corr_per_row, 0.8)

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
            .replace(".txt", ""),  # Extract the part after "HLA_" and before ".txt"
        )
        
        # print(highest_corr_per_row_sorted)
        for row, (col, corr) in highest_corr_per_row_sorted:
            # Generate title based on row and column
            title = f"{row} -> {col}: {corr:.2f}"
            # title2 = f"naturally presented logo of{self.formate_HLA_DB(col)} -> {row}: {corr:.2f}"
            title2 = f"naturally presented logo of {self.formate_HLA_DB(col)} -> {row}: {corr:.2f}"

            # Find the image paths
            image_path = self._find_gibbs_image_path(row.split(".")[1], image_folder)
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
            self.console.log("Error: No images were found to create a grid.", style="bold red")
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

    ## New Image Grid Module
    def generate_image_grid(corr):
        """
        Generate an image grid for the correlation results.
        """
        pass


def run_cluster_search(args):
    # if args.output_folder is None:
    #     os.makedirs("output", exist_ok=True)
    #     output_folder = "output"
    cluster_search = ClusterSearch()
    cluster_search.console.rule(
        "[bold red]Stage 1/4: Search all cluster best match HLA type."
    )

    cluster_search.compute_correlations(
        args.gibbs_folder,
        args.reference_folder,
        args.n_clusters,
        args.output,
        args.hla_types,
    )

    cluster_search.console.rule(
        "[bold red]Stage 2/4: Finding best matching Naturally presented HLA ."
    )
    # if args.output_folder is None:
    cluster_search.plot_heatmap(args.output)

    cluster_search.console.rule("[bold red]Stage 3/4: Cheking the HLA.")

    cluster_search.check_HLA_DB(args.hla_types, args.reference_folder)
    # cluster_search.create_image_grid(cluster_search.correlation_dict, os.path.join(args.gibbs_folder, 'logos'), os.path.join(args.reference_folder, 'images'), os.path.join(args.output_folder, 'image_grid_D90.png'), HLA_list=cluster_search.valid_HLA)

    cluster_search.console.rule(
        "[bold red]Stage 4/4: Finding Best Matched HLA for Each cluster."
    )
    cluster_search.create_image_grid(
        cluster_search.correlation_dict,
        os.path.join(args.gibbs_folder, "logos"),
        str(args.reference_folder).replace("/output_matrices_human", "").replace("output_matrices",""),
        os.path.join(args.output, "compare_motif_corr.png"),
        HLA_list=cluster_search.valid_HLA,
    )
    # logging.info("Process completed successfully.")
    if args.log:
        log_file_path = os.path.join(cluster_search._outfolder, "search_cluster.log")
        save_console_log()
        # raise FileNotFoundError(f"Log file not found: {log_file_path}")

## Remove thie after test

if __name__ == "__main__":
    # ClusterSearch().compute_correlations(
    #     "data/9mersonly",
    #     "data/ref_data/Gibbs_motifs_mouse/output_matrices",
    #     "all",
    #     "data/outputM",
    # )
        
    # ClusterSearch()._compute_and_log_correlation(
    #     "data/9mersonly",
    #     "data/ref_data/Gibbs_motifs_human/output_matrices_human",
    #     "cluster_1of5.mat",
    #     "HLA_A_02_01.txt",
    # )
    # print(sys.argv)
    
    print(ClusterSearch()._db_loader("data/ref_data/","mouse"))
        
    # run_cluster_search(
    #     argparse.Namespace(
    #         credits=False,
    #         gibbs_folder="data/9mersonly",
    #         hla_types=None,
    #         log=False,
    #         n_clusters="5",
    #         output="data/outputM",
    #         processes=4,
    #         reference_folder="data/ref_data/Gibbs_motifs_mouse/output_matrices",
    #         version=False,
    #     )
    # )