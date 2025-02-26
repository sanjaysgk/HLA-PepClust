"""Error logs for the CLI."""


def _file_not_exists(file_path: str):
    """Error message for file not found."""
    return f"File not found: {file_path}"


def _png_not_exists(file_path: str):
    """Error message for PNG file not found."""
    return f"PNG file not found: {file_path}"
