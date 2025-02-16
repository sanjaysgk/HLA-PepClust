"""HLA-PEPCLUST: CLI Application for HLA Motif Discovery"""

import argparse
import sys
from multiprocessing import Pool
from rich.text import Text
from rich_argparse import RichHelpFormatter
from . import __version__
from cli.cluster_search import run_cluster_search
from cli.logger import CONSOLE

# from pyfiglet import Figlet
# fig = Figlet(font="slant")
# tool_icon= fig.renderText("HLA-PEPCLUST")
# print(ascii_text)


def _welcome():
    """Display application banner."""

    tool_icon = """
    __  ____    ___         ____  __________  ________    __  _____________
   / / / / /   /   |       / __ \/ ____/ __ \/ ____/ /   / / / / ___/_  __/
  / /_/ / /   / /| |______/ /_/ / __/ / /_/ / /   / /   / / / /\__ \ / /   
 / __  / /___/ ___ /_____/ ____/ /___/ ____/ /___/ /___/ /_/ /___/ // /    
/_/ /_/_____/_/  |_|    /_/   /_____/_/    \____/_____/\____//____//_/     
                                                                           
    """

    CONSOLE.print(tool_icon, style="blue")


def _print_credits(credits=False):
    """Print software credits to terminal."""
    text = Text()
    text.append("\n")
    if credits:
        text.append("Please cite: \n", style="bold")
        text.append(
            "GibbsCluster - 2.0 (Simultaneous alignment and clustering of peptide data)\n",
            style="bold link https://services.healthtech.dtu.dk/services/GibbsCluster-2.0/",
        )
        text.append(
            "Seq2Logo - 2.0 (Visualization of amino acid binding motifs)\n",
            style="bold link https://services.healthtech.dtu.dk/services/Seq2Logo-2.0/",
        )
        text.append(
            "MHC Motif Atlas (MHC motif PSM matrics are genrated using mhcmotifatlas please cite [link])\n\n",
            style="bold link http://mhcmotifatlas.org/home",
        )
    else:
        text.append(
            "HLA-PEPCLUST", style="bold link https://www.monash.edu/research/compomics/"
        )
        text.append(f" (v{__version__})\n", style="bold")
    if credits:
        text.append(
            "HLA motif finder pipeline for identifying peptide motif immunopeptididomic data.\n",
            style="italic",
        )
    text.append("Developed at Li Lab / Purcell Lab, Monash University, Australia.\n")
    text.append("Please cite: ")
    if credits:
        text.append(
            "Sanjay Krishna, Nathon Craft & Chen Li et al. bioRxiv (2024)",
            style="link https://www.monash.edu/research/compomics/",
        )
    else:
        text.append(
            "Sanjay Krishna & Chen Li et al. bioRxiv (2024)",
            style="ttps://www.monash.edu/research/compomics/",
        )
    text.append("\n")
    if credits:
        text.stylize("#006cb5")
    CONSOLE.print(text)


def main():
    """Main function to parse CLI arguments and execute the pipeline."""
    _welcome()
    _print_credits()
    parser = argparse.ArgumentParser(
        description="HLA-PEPCLUST: Identify peptide motifs in immunopeptidomics data 🧬🧬🧬.",
        formatter_class=lambda prog: RichHelpFormatter(prog, max_help_position=42),
    )

    parser.add_argument(
        "gibbs_folder", type=str, help="Path to the test folder containing matrices."
    )
    parser.add_argument(
        "reference_folder",
        type=str,
        help="Path to the human reference folder containing matrices.",
    )

    parser.add_argument(
        "-hla",
        "--hla_types",
        nargs="*",
        default=None,
        help="List of HLA types to search (defaults to all).",
    )
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel processes to use.",
    )
    parser.add_argument(
        "-n",
        "--n_clusters",
        type=str,
        default="all",
        help="Number of clusters to search for.",
    )
    parser.add_argument(
        "-k", "--best_KL", action="store_true", help="Find the best KL divergence only."
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output", help="Path to the output folder."
    )
    parser.add_argument("-l", "--log", action="store_true", help="Enable logging.")

    parser.add_argument(
        "-c",
        "--credits",
        action="store_true",
        help="Show credits for the motif database pipeline.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s [blue]version[/blue] [i][blue]{__version__}[/blue][/i]",
        help="Show the version of the pipeline.",
    )

    args = parser.parse_args()

    if args.credits:
        _print_credits(credits=True)
        sys.exit(0)

    if args.processes > 1:
        with Pool(args.processes) as pool:
            pool.map(run_cluster_search, [args])
    else:
        run_cluster_search(args)


if __name__ == "__main__":
    main()
