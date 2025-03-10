"""
Logging configuration for HLA-PEPCLUST.

This module provides centralized logging for the application with
rich formatting and optional file output capabilities.
"""

import logging
import os
from rich.logging import RichHandler
from rich.console import Console
from rich.table import Table
import click

# Console for Rich library
CONSOLE = Console(record=True)

# Logging level mapping
LOG_LEVELS = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


def configure_logging(level: str = "info", log_to_file: bool = False, log_file: str = None) -> logging.Logger:
    """
    Configures logging with RichHandler and optional file logging.

    Args:
        level (str): Logging level ("critical", "error", "warning", "info", "debug"). Defaults to "info".
        log_to_file (bool): If True, logs are also saved to a file. Defaults to False.
        log_file (str, optional): Path to log file. If None, defaults to "cluster_search.log" in current directory.

    Returns:
        logging.Logger: Configured logger instance.
    """
    log_level = LOG_LEVELS.get(level.lower(), logging.INFO)
    logger = logging.getLogger("hla_pepclust")

    # Configure handlers
    handlers = [
        RichHandler(
            rich_tracebacks=True,
            tracebacks_show_locals=True,
            tracebacks_suppress=[click],
        )
    ]
    
    if log_to_file:
        if log_file is None:
            log_file = "cluster_search.log"
        
        # Ensure the directory exists
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir)
            
        file_handler = logging.FileHandler(log_file)
        handlers.append(file_handler)

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="[%X]",
        handlers=handlers,
    )

    logger.setLevel(log_level)
    return logger


def save_console_log(file_name: str = "cluster_search.log"):
    """
    Saves console logs to a specified file.

    Args:
        file_name (str): The filename to save the log. Defaults to "cluster_search.log".
    """
    log_dir = os.path.dirname(file_name)
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
        
    with open(file_name, "w") as log_file:
        log_file.write(CONSOLE.export_text())


def display_search_results(highest_co_mat: dict, threshold: float):
    """
    Displays search results in a formatted Rich table.

    Args:
        highest_co_mat (dict): Dictionary containing search results.
        threshold (float): Score threshold for highlighting results.
    """
    table = Table(title="Search Results")
    table.add_column("Cluster", style="magenta")
    table.add_column("Best HLA Match", justify="center")
    table.add_column("Score", justify="right")

    for key, value in highest_co_mat.items():
        score = value[1]
        score_color = (
            "green" if score >= threshold
            else "yellow" if score >= 0.5 and score < threshold else "red"
        )
        best_hla_color = "green" if score >= threshold else "gray"
        table.add_row(
            f"[{best_hla_color}]{str(key).split('/')[-1]}[/]",
            f"[{best_hla_color}]{str(value[0]).split('/')[-1].replace('.txt', '')}[/]",
            f"[{score_color}]{score:.3f}[/]",
        )

    CONSOLE.print(table)