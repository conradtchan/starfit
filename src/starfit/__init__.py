import os
from pathlib import Path

DATA_DIR = Path(__file__).parent.resolve() / "data"
del Path

# If user has defined an environment variable for custom data files, use that instead
user_data_dir = os.getenv("STARFIT_DATA")
if user_data_dir is not None:
    DATA_DIR = user_data_dir

from importlib import metadata

__version__ = metadata.version("starfit")


from .fit import Direct, Double, Single
from .ga import Ga

__all__ = ["Direct", "Double", "Single", "Ga"]
