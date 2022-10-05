"""
provide specific abundance sets
"""

import os
import os.path
import re
import time

import numpy as np

from .. import DATA_DIR
from . import isotope
from .abumodel import AbuModel
from .abuset import AbuSet
from .utils import CachedAttribute


class SolAbu(AbuSet):
    """
    Special abundance set with interface to load from disk files.

    Most of the functionallity is in AbuSet, though.

    """

    default_sets = {
        "AG89": "solag89.dat",
        "GN93": "solgn93.dat",
        "Lo03": "sollo03.dat",
        "Lo09": "sollo09.dat",
        "As09": "solas09.dat",
        "As12": "solas12.dat",
        "As09s": "solas09_sol_surf_present.dat",
    }

    default = "As09s"

    def __init__(self, name=None, silent=False):
        """
        Create abundance from set name.
        """
        super().__init__(silent=silent)

        if name is None:
            name = self.default

        self.iso = np.array([], dtype=object)
        self.abu = np.array([], dtype=np.float64)
        self.comment = ()
        self._load_abu(self.path(name, self.default_sets))
        self.normalize()
        self.is_sorted = True
        self.sort()

    def path(self, name, sets):
        """
        Find file from path or set name.

        TODO - should probably be named differently
        """
        self.setup_logger(silent=self.silent)
        if name in sets:
            self.name = name
            path = os.path.join(DATA_DIR, "ref")
            filename = os.path.join(path, sets[name])
        else:  # HERE we should add error treatment
            self.name = name
            filename = name
        self.close_logger()
        return filename

    def _load_abu(self, filename):
        self._from_dat(filename, silent=self.silent)


class BBNAbu(SolAbu):
    """
    load BBN set
    """

    default_sets = {"F02": "bbnf02.dat"}
    default = "F02"


class ScaledSolar(AbuModel):
    """
    Special abundance set created from scaled solar abundance set.

    """

    version = "10000"

    def _abu_massfrac_raw(self, scale):
        """
        Raw scaled solar abundances
        """

        scaled = self.sun * scale + self.bbn * (1 - scale)

        # beyond-solar scaling
        if scale > 1.0:
            (jj,) = np.argwhere(scaled.iso == isotope.Ion("He4"))
            # make sure we have same set of istopes
            self.bbn = (self.sun * 0) + self.bbn
            for j in np.argwhere(scaled.abu < self.sun.abu).flat:
                scaled.abu[jj] += scaled.abu[j]
                scaled.abu[j] = self.sun.abu[j] * np.exp(
                    (scale - 1) * (1 - self.bbn.abu[j] / self.sun.abu[j])
                )
                scaled.abu[jj] -= scaled.abu[j]
        scaled.normalize()

        return scaled.abu

    # use same defaults and definitions as in SolAbu
    def __init__(self, scale=1, **kw):

        """
        Create abundance from set name.

        Use simple algorithm:
        X = X_sun * scale + X_BBN * (1 - scale)

        For stuff that get less, assume exponential decrease beyond
        solar abundance.  The idea is that less can get incorporated
        in stars.  The difference goes into He since H1 itself
        decreases as well.
        """

        silent = kw.setdefault("silent", False)
        self.setup_logger(silent=silent)

        solar = kw.pop("solar", None)
        if solar is None:
            solar = SolAbu.default
        zero = kw.pop("zero", None)
        if zero is None:
            zero = BBNAbu.default

        self.sun = SolAbu(solar)
        self.bbn = BBNAbu(zero)

        check = kw.get("check", False)
        if check:
            assert len(self.bbn) == len(self.sun) == len(self.ions)

        super().__init__(scale, **kw)

        self.is_sorted = self.sun.is_sorted
        self.comment = (
            f"Version {self.version:6s} - {time.asctime(time.gmtime()) + ' UTC':s}",
            f"Scaled solar abundances: {scale:g} solar",
            f"Sun: {solar:s} - {self.sun.filename:s}",
            f"BBN: {zero:s} - {self.bbn.filename:s}",
            "X = {:8G}, Y = {:8G}, Z = {:8G}".format(*self.XYZ()),
        )
        self.logger.info("X = {:8G}, Y = {:8G}, Z = {:8G}".format(*self.XYZ()))
        self.close_logger()


class ScaledSolarHelium(ScaledSolar):
    """
    Special abundance set created from scaled solar abundance set and
    overwrite He abundance.

    Version history:
    10000: created
    10001: add light isotope scaling
    """

    version = "10001"

    def _abu_massfrac_raw(self, scale):
        """
        Raw scaled solar abundances.
        """
        # TODO - use abuset object for abu

        abu = super()._abu_massfrac_raw(scale)

        if self.helium is None:
            return abu

        jH1, jHe4 = self.sun.index(("H1", "He4"))

        if not self.scale_light:
            abu[jH1] -= self.helium - abu[jHe4]
            abu[jHe4] = self.helium
        else:
            # do general scaling of light isotopes.
            # most primitive: just assume the are destroyed
            # along with H1 (anything else is model extrapolation)
            # for now we only do H2, He3 as extra light isotopes
            #
            # Assume the following reactions:
            #    2 H2  --> He4
            #    2 He3 --> He4 + 2 H1

            jH2, jHe3 = self.sun.index(("H2", "He3"))
            # jLi6 = np.argwhere(self.sun.iso == isotope.Ion('Li6'))
            # jLi7 = np.argwhere(self.sun.iso == isotope.Ion('Li7'))
            # jBe9 = np.argwhere(self.sun.iso == isotope.Ion('Be9'))
            # jB10 = np.argwhere(self.sun.iso == isotope.Ion('B10'))
            # jB11 = np.argwhere(self.sun.iso == isotope.Ion('B11'))

            xH1 = abu[jH1]
            xH2 = abu[jH2]
            xHe3 = abu[jHe3]
            xHe4 = abu[jHe4]

            dHe4 = self.helium - xHe4
            source = xH1 + xH2 + 2 / 3 * xHe3
            f = dHe4 / source

            abu[jH1] -= f * (xH1 - xHe3 / 3)
            abu[jH2] -= f * xH2
            abu[jHe3] -= f * xHe3
            abu[jHe4] = self.helium

        return abu

    # use same defaults and definitions as in SolAbu
    def __init__(self, scale=1, helium=None, scale_light=False, **kwargs):

        (
            """
        Based on scaled solar, overwrite He4, rest in H1.

        parameters:
          scale:
            scale on solar abundance
          helium:
            desired He4 mass fraction
          scale_light:
            also scale other light istopes asumed to be
            destroyed to same fraction as H1

        Keeps H2 and He3 values by default.

        Base class documnetation:
        -------------------------
        """
            + super().__doc__
        )
        self.helium = helium
        self.scale_light = scale_light
        super().__init__(scale, **kwargs)
        if helium is None:
            helium = self("He4")
        self.comment += (f"Helium set to mass fraction: : {helium:g}",)


class Asplund2009Data(AbuSet):
    """
    Routine to load Asplund 2009 solar abundances

    Here is what I used to call:
    x = Asplund2009Data()
    x.write_dat('~/kepler/local_data/solas12.dat')
    x.write_bg('~/kepler/local_data/solas12g')
    """

    def __init__(
        self,
        filename="~/Plots/solar/Asplund2009-isotopes_protosun.dat",
        comment=None,
        silent=False,
    ):
        """
        Load abundace set from Aspund "dat" file.

        TODO - add option to show comment
        """
        self.setup_logger(silent=silent)
        comment = tuple(comment)

        xre = re.compile("[-+a-zA-Z0-9.]+")
        iso = np.array([], dtype=object)
        abu = np.array([], dtype=np.float64)

        filename = os.path.expanduser(filename)

        with open(filename, "r") as f:
            self.logger_file_info(f)
            comment += (
                "",
                f'Generated from file "{filename:s}".',
                "Original file comments follow:",
                "",
            )
            for line in f:
                if not line.startswith((";", "#")):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum == 0:
                        continue
                    if xnum == 5:
                        xiso = isotope.Ion(A=int(xdata[2]), Z=int(xdata[1]))
                        xabu = np.double(xdata[4])
                    else:
                        print(line)
                        raise IOError("bad format")
                    iso = np.append(iso, xiso)
                    abu = np.append(abu, xabu)
                else:
                    comment += (line[2:].rstrip(),)
        m = Mass()

        # well, this could require tests...
        abu = np.array([a * m(i) for i, a in zip(iso, abu)])
        abu = abu / abu.sum()

        super().__init__(iso=iso, abu=abu, comment=comment)
        message = f"{iso.size:3d} isotopes loaded in"
        self.close_logger(timing=message)


from .logged import Logged


# seeems to be not used at present
# add option to AbuSet etc. to use proper masses
class Mass(Logged):
    """
    Object to hold ion masses.

    Not clear at this point where the informastion will come from.
    Currently load audi 2003

    For the rest: for now just use A

    TODO: use N*m_n + Z*m_p - Q/c**2 (from bdat)
    TODO: Option to just use A
    """

    default_data = "Au03"

    def __init__(self, data=default_data, silent=False):
        """
        TODO - implement data
        """
        self.setup_logger(silent=silent)
        path = os.getenv("KEPLER_DATA")
        if not path:
            path = os.path.join(os.path.expanduser("~"), "kepler", "local_data")
            self.logger.warning("using default path " + path)
        filename = os.path.join(path, "masses_audi_2003.dat")

        self.comment = ()
        self.iso = np.array([], dtype=object)
        self.mass = np.array([], dtype=np.float64)

        xre = re.compile("[-+a-zA-Z0-9.]+")
        with open(filename, "r") as f:
            self.logger_file_info(f)
            for line in f:
                if not line.startswith((";", "#")):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum == 0:
                        continue
                    if xnum == 2:
                        xion, xabu = tuple(xdata)
                    else:
                        print(line)
                        raise IOError("bad format")
                    self._append(isotope.Ion(xion), np.double(xabu))
                else:
                    self.comment += (line[2:],)
        message = f"{len(self.iso):3d} masses loaded in"
        self.close_logger(timing=message)

    def _append(self, iso, abu):
        self.iso = np.append(self.iso, iso)
        self.mass = np.append(self.mass, abu)

    def append(self, iso, abu):
        self._append(iso, abu)
        del self.A
        del self.Z
        del self.N
        del self.DM
        del self.BE

    def __getitem__(self, ion):
        """
        Return mass of ion in amu (or approximate)

        TODO: accept list
        """
        try:
            (i,) = np.argwhere(self.iso == ion)
            return self.mass[i[0]]
        except:
            pass
        return np.double(isotope.Ion(ion).A)

    def __str__(self):
        return (
            "mass("
            + ", ".join(
                [f"{iso.Name():s}: {mass:f}" for iso, mass in zip(self.iso, self.mass)]
            )
            + ")"
        )

    __repr__ = __str__

    @CachedAttribute
    def A(self):
        return isotope.Ion.ufunc_A(self.iso)

    @CachedAttribute
    def Z(self):
        return isotope.Ion.ufunc_Z(self.iso)

    @CachedAttribute
    def N(self):
        return isotope.Ion.ufunc_N(self.iso)

    @CachedAttribute
    def DM(self):
        """mass excess in amu (not correct unit)"""
        return self.mass - self.A

    @CachedAttribute
    def BE(self):
        """binding energy in amu (not correct unit)"""
        mn = 1.008664916
        mp = 1.0078250321
        return self.mass - self.Z * mp - self.N * mn
