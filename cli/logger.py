import logging
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


def configure_logging(level: str = "info", log_to_file: bool = False) -> logging.Logger:
    """
    Configures logging with RichHandler and optional file logging.

    Args:
        level (str): Logging level ("critical", "error", "warning", "info", "debug"). Defaults to "info".
        log_to_file (bool): If True, logs are also saved to a file. Defaults to False.

    Returns:
        logging.Logger: Configured logger instance.
    """
    log_level = LOG_LEVELS.get(level.lower(), logging.INFO)
    logger = logging.getLogger(__name__)

    # Configure handlers
    handlers = [
        RichHandler(
            rich_tracebacks=True,
            tracebacks_show_locals=True,
            tracebacks_suppress=[click],
        )
    ]
    if log_to_file:
        file_handler = logging.FileHandler("cluster_search_pipeline.log")
        handlers.append(file_handler)

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
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
    with open(file_name, "w") as log_file:
        log_file.write(CONSOLE.export_text())


def display_search_results(highest_co_mat: dict, threshold: int):
    """
    Displays search results in a formatted Rich table.

    Args:
        highest_co_mat (dict): Dictionary containing search results.
        threshold (int): Score threshold for highlighting results.
    """
    table = Table(title="Search Results")
    table.add_column("Cluster", style="magenta")
    table.add_column("Best HLA Match", justify="center")
    table.add_column("Score", justify="right")

    for key, value in highest_co_mat.items():
        score = value[1]
        score_color = (
            "green"
            if score >= threshold
            else "yellow" if score > 0.5 < threshold else "red"
        )
        best_hla_color = "green" if score >= threshold else "gray"
        table.add_row(
            f"[{best_hla_color}]{str(key).split('/')[-1]}[/]",
            f"[{best_hla_color}]{str(value[0]).split('/')[-1].replace('.txt', '')}[/]",
            f"[{score_color}]{score}[/]",
        )

    CONSOLE.print(table)
