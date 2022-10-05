"""
utilities for kepler tools
"""

from .utils import iterable


def mass_equal(mass1, mass2):
    """
    Compare mass values within round-off.
    (strings are converted to float)
    """
    if isinstance(mass1, str):
        mass1 = float(mass1)
    if isinstance(mass2, str):
        mass2 = float(mass2)
    return abs(mass1 - mass2) / (mass1 + mass2) < 1.0e-12


def mass_string(masses, decimals=None, powers=False):
    """
    convert mass number to string
    """
    masses = iterable(masses)
    if decimals is None:
        decimals = 0
    xmass = []
    for mass in masses:
        if isinstance(mass, str):
            mass = float(mass)
        xm = f"{mass:.6f}".rstrip("0").rstrip(".")
        if decimals > 0:
            if xm.count(".") == 0:
                xm += "." + "0" * decimals
        if powers:
            if xm.count(".") == 0:
                if xm.endswith("0" * 24):
                    xm = xm[:-24] + "Y"
                elif xm.endswith("0" * 21):
                    xm = xm[:-21] + "Z"
                elif xm.endswith("0" * 18):
                    xm = xm[:-18] + "E"
                elif xm.endswith("0" * 15):
                    xm = xm[:-15] + "P"
                elif xm.endswith("0" * 12):
                    xm = xm[:-12] + "T"
                elif xm.endswith("0" * 9):
                    xm = xm[:-9] + "G"
                elif xm.endswith("0" * 6):
                    xm = xm[:-6] + "M"
                elif xm.endswith("0" * 3):
                    xm = xm[:-3] + "k"
        xmass += [xm]
    if len(xmass) == 1:
        xmass = xmass[0]
    return xmass


def mass_formatter(*args, **kwargs):
    """
    function to format mass string for use with ticker
    """
    return mass_string(args[0])


import matplotlib.ticker


class MassFormatter(matplotlib.ticker.FuncFormatter):
    def __init__(self, *args, **kwargs):
        super().__init__(mass_formatter)


class MissingModels(Exception):
    """
    Exception raised for KEPLER data files missing models in sequence
    """

    def __init__(self, models, filename):
        self.models = models
        self.filename = filename

    def __str__(self):
        return f"Missing models in file {self.filename}: " + ", ".join(
            [f"{x:d}" for x in self.models]
        )


class RecordVersionMismatch(Exception):
    """
    Exception raised for KEPLER data files with different versions.
    """

    def __init__(self, models, filename):
        self.models = models
        self.filename = filename

    def __str__(self):
        return f"Missing models in file {self.filename}: " + ", ".join(
            [f"{x:d}" for x in self.models]
        )


class UnkownVersion(Exception):
    def __init__(self, version=None, record=None):
        self.version = version
        self.record = record
