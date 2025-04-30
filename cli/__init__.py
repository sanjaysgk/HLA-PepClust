__version__ = "1.1.0-beta"


from warnings import filterwarnings
from rich.traceback import install

install(show_locals=True)  # type: ignore

# Suppress warnings for missing database path
filterwarnings(
    "ignore",
    message="Database path is not set. Using default path",
    category=UserWarning,
    module="cli.database_gen",
    append=True
)
