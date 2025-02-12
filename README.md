# HLA-PepClust

## Introduction

HLA-PepClust `CLI` module is designed for identifying HLA type by clustering peptide sequences based on their
## Prerequisites

Ensure you have the following installed on your system:
- Python 3.9 or higher
- `pip` (Python package installer)
## Download or clone the git repo
    ```bash
    git clone https://github.com/Sanpme66/HLA-PepClust.git
    ```

## Setting Up the Python Environment

1. **Create a virtual environment**:
    ```bash
    python3 -m venv hlapepclust-env
    ```

2. **Activate the virtual environment**:
    - On macOS and Linux:
      ```bash
      source hlapepclust-env/bin/activate
      ```
    - On Windows:
      ```bash
      .\hlapepclust-env\Scripts\activate
      ```

3. **Upgrade `pip`**:
    ```bash
    pip install --upgrade pip
    ```

## Installing Dependencies

1. **Navigate to the project directory**:
    ```bash
    cd HLA-PepClust/
    ```

2. **Install the required packages**:
    ```bash
    pip install -r requirements.txt
    ```

## Running HLA-PepClust

1. **Run the main script**:

    ```bash
    python CLI/cluster_search.py <input_data_path> <reference_data_path> --hla_types <hla_types> --n_clusters <number_of_clusters> --output <output_path>
    ```
    ### Command Line Arguments

    The `cluster_search.py` script accepts the following command line arguments:

    - `gibbs_folder` (str): Path to the test folder containing matrices.
    - `reference_folder` (str): Path to the human reference folder containing matrices.
    - `--output_folder` (str, optional): Path to the output folder. Default is "output".
    - `--hla_types` (list, optional): List of HLA types to search (defaults to all if none specified).
    - `--processes` (int, optional): Number of processes to use. Default is 4.
    - `--n_clusters` (str, optional): Number of clusters to search for. Default is "all".
    - `--best_KL` (bool, optional): Find the best KL divergence only. Default is False.
    - `--output` (str, optional): Path to the output folder. Default is "output".

    Example usage:
    ```bash
    python CLI/cluster_search.py data/D90_HLA_3844874 data/ref_data/Gibbs_motifs_human/output_matrices_human --hla_types A0201,A0101,B1302,B3503,C0401 --n_clusters 6 --output test_results --processes 4
    ```
    

## Deactivating the Virtual Environment

After you are done, you can deactivate the virtual environment by running:
```bash
deactivate
```

## Conclusion

More detailed instruction comming soon......!!!!!!
