# HLA-PEPCLUST

HLA Motif Finder for Immunopeptidomics Data

## Overview

HLA-PEPCLUST is a command-line tool designed to identify human leukocyte antigen (HLA) or major histocompatibility complex (MHC) binding motifs in immunopeptidomic data. The pipeline compares clustering results from GibbsCluster with a database of known HLA/MHC binding motifs to identify the best matches.

## Installation

```bash
# Clone the repository
git clone https://github.com/Sanpme66/HLA-PepClust.git
cd HLA-PepClust

# Create and activate a virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Database Generation

Before running the cluster search, you need to generate the reference database:

```bash
python -m cli.main -db path/to/config.json
```

### Cluster Search

To run the cluster search, use the following command:

```bash
python -m cli.main path/to/gibbs_results path/to/reference_folder -o output/path -s species -hla HLA1,HLA2 -n num_clusters
```
```bash
python -m cli.main data/D105heartW6:32_361712_originally_labeled_D120/ data/ref_data/Gibbs_motifs_human -o output/resultsD105 -s human -n 6
```

#### Parameters:

- `gibbs_folder`: Path to the folder containing GibbsCluster results
- `reference_folder`: Path to the reference database folder
- `-o, --output`: Path to save output (default: "output")
- `-s, --species`: Species to search (human or mouse, default: human)
- `-hla, --hla_types`: List of HLA/MHC allotypes to search (comma-separated, defaults to all)
- `-n, --n_clusters`: Number of clusters to search for (default: "all")
- `-t, --threshold`: Threshold for motif similarity (default: 0.5)
- `-p, --processes`: Number of parallel processes to use (default: 4)
- `-im, --immunolyser`: Enable immunolyser output
- `-l, --log`: Enable logging
- `-c, --credits`: Show detailed credits and citations
- `-v, --version`: Show the version of the pipeline

### Example

```bash
python -m cli.main data/gibbs_output data/ref_data -o output/results -p 6 -hla A0101,A0201,B0702 -s human -n 6 --log -im
```

## Output

HLA-PEPCLUST generates the following outputs:

- HTML report visualizing the correlation between clusters and reference motifs
- Heatmap of correlation values
- Interactive plots of amino acid correlations
- Tables of matching scores
- Log files (if logging is enabled)

## Citation

If you use HLA-PEPCLUST in your research, please cite:

Sanjay Krishna, Nathon Craft & Chen Li et al. bioRxiv (2024)

## Dependencies

- pandas
- numpy
- matplotlib
- seaborn
- altair
- Pillow
- jinja2
- rich

## License

MIT License

## Contact

For questions or issues, please contact:
- Sanjay Krishna (sanjay.krishna@monash.edu)
- Chen Li (chen.li@monash.edu)
- Prithvi