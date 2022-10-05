from pathlib import Path

DATA_DIR = Path(__file__).parent.resolve() / "data"
del Path

from .fit import Direct, Double, Single
from .ga import Ga

__all__ = ["Direct", "Double", "Single", "Ga"]
