from os import getenv
from pathlib import Path

DATA_DIR = Path(__file__).parent.resolve() / "data"

# If user has defined an environment variable for custom data files, use that instead
user_data_dir = getenv("STARFIT_DATA")
if user_data_dir is not None:
    DATA_DIR = Path(user_data_dir).exanduser().resolve()

SOLAR = "sollo22.dat"
user_solar = getenv("STARFIT_SOLAR")
if user_solar is not None:
    SOLAR = user_solar

from importlib import metadata

__version__ = metadata.version("starfit")


from .fit import Direct, Double, Single
from .ga import Ga

__all__ = ["Direct", "Double", "Single", "Ga"]

del Path
