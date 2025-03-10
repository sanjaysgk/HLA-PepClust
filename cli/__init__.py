# This file is part of the HLA-PEPCLUST software package.
from warnings import filterwarnings
from rich.traceback import install

install(show_locals=True)  # type: ignore

"""
HLA-PEPCLUST - HLA Motif Finder for Immunopeptidomic Data

A Python package for identifying peptide motifs in immunopeptidomic data 
and matching them to known HLA/MHC allotypes.
"""

__version__ = "1.0.0-dev"
__author__ = "Sanjay Krishna,Prithvi Munday,Chen Li"
__email__ = "sanjay.sondekoppagopalakrishna@monash.edu"
