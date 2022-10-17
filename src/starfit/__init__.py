import os
from pathlib import Path

DATA_DIR = Path(__file__).parent.expanduser().resolve() / "data"
DATA_DIRS = [DATA_DIR]

# If user has defined an environment variable for custom data files,
# add that to the list of data directories to search, assuming it is valid
user_data_dir = os.getenv("STARFIT_DATA")
if user_data_dir:
    user_data_dir = Path(user_data_dir).expanduser().resolve()
    if user_data_dir.is_dir():
        DATA_DIRS = [user_data_dir, DATA_DIR]
    else:
        print(f' [StarFit] STARFIT_DATA="{user_data_dir}" is not valid.  Ignoring.')

SOLAR = "sollo22.dat"
user_solar = os.getenv("STARFIT_SOLAR")
if user_solar is not None:
    SOLAR = user_solar

BBN = "bbnc19.dat"
user_bbn = os.getenv("STARFIT_BBN")
if user_bbn is not None:
    BBN = user_bbn

from importlib import metadata

__version__ = metadata.version("starfit")


from .fit import Direct, Double, Single
from .ga import Ga

__all__ = ["Direct", "Double", "Single", "Ga"]

del Path, os
del user_data_dir, user_solar, user_bbn
