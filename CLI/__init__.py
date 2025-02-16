__version__ = "0.0.1-dev"


from warnings import filterwarnings
from rich.traceback import install
install(show_locals=True)  # type: ignore

# mzmlb is not used, so hdf5plugin is not needed
filterwarnings(
    "ignore",
    message="hdf5plugin is missing",
    category=UserWarning,
    module="psims.mzmlb",
)
