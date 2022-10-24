from os import getenv
from pathlib import Path

DATA = "data"
DATA_DIR = Path(__file__).parent.expanduser().resolve() / DATA
DATA_DIRS = [DATA_DIR]

DB = "db"
REF = "ref"
STARS = "stars"

# If user has defined an environment variable for custom data files,
# add that to the list of data directories to search, assuming it is valid
user_data_dir = getenv("STARFIT_DATA")
if user_data_dir:
    for user_data_dir in user_data_dir.split(":"):
        user_data_dir = Path(user_data_dir).expanduser().resolve()
        if user_data_dir.is_dir():
            DATA_DIRS = [user_data_dir] + DATA_DIRS
        else:
            print(f' [StarFit] STARFIT_DATA="{user_data_dir}" is not valid.  Ignoring.')

SOLAR = "sollo22.dat"
user_solar = getenv("STARFIT_SOLAR")
if user_solar is not None:
    SOLAR = user_solar

BBN = "bbnc19.dat"
user_bbn = getenv("STARFIT_BBN")
if user_bbn is not None:
    BBN = user_bbn

from importlib import metadata

__version__ = metadata.version("starfit")


from .direct import Direct
from .ga import Ga
from .multi import Multi
from .single import Single
from .star import Star

__all__ = ["Single", "Multi", "Ga", "Direct", "Star"]

del Path, getenv
del user_data_dir, user_solar, user_bbn
