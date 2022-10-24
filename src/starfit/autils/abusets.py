"""
provide specific abundance sets
"""

# standard packages
import os
import os.path
import re
import sys
import time
from pathlib import Path
from urllib import request
from urllib.error import HTTPError

# non-standard packages
import numpy as np

# alex's packages
from .abumodel import AbuModel
from .abuset import AbuSet
from .human import byte2human
from .isotope import ion as I
from .isotope import ufunc_A, ufunc_E, ufunc_N, ufunc_Z
from .logged import Logged
from .utils import CachedAttribute, stuple

# path = Path(__file__).parent.resolve() / ".abusets"
path = Path(__file__).parent.resolve() / "../data/ref"
page = "https://2sn.org/Download/abusets/"


def downloadfile(filename):
    url = page + filename
    try:
        response = request.urlopen(url)
    except HTTPError as e:
        print(f"[abusets] {e}: {url}")
        return
    filename = Path(filename).expanduser()
    if not filename.is_absolute():
        filename = path / filename
    print(f" [abusets.downloadfile] downloading {url}")
    content = response.read()
    print(f" [abusets.downloadfile] writing {byte2human(len(content))} to {filename!s}")
    filename.parent.mkdir(parents=True, exist_ok=True)
    filename.write_bytes(content)


# these files can be replaced by user in the .abusets directory
_abusets = (
    "solag89.dat",
    "solgn93.dat",
    "sollo03.dat",
    "sollo09.dat",
    "solas09.dat",
    "solas12.dat",
    "sollo20.dat",
    "sollo22.dat",
    "bbnf02.dat",
    "bbnc16.dat",
    "bbnc19.dat",
)

for abuset in _abusets:
    if not (path / abuset).exists():
        downloadfile(abuset)


class SolAbu(AbuSet):
    """
    Special abundance set with interface to load from disk files.

    Most of the functionality is in AbuSet, though.

    """

    default_sets = {
        "AG89": "solag89.dat",
        "GN93": "solgn93.dat",
        "Lo03": "sollo03.dat",
        "Lo09": "sollo09.dat",
        "As09": "solas09.dat",  # Asplund 09, current sun
        "As12": "solas12.dat",  # Asplund 09, presolar (pc)
        "Lo20": "sollo20.dat",  # Lodders 20, presolar
        "Lo22": "sollo22.dat",  # Lodders 20, current sun (reconstructed by Alex)
    }

    default = "Lo20"

    def __init__(
        self,
        name=None,
        silent=False,
        **kwargs,
    ):
        """
        Create abundance from set name.
        """
        super().__init__(silent=silent)

        if name is None:
            name = self.default
        self.iso = np.array([], dtype=object)
        self.abu = np.array([], dtype=np.float64)
        self.comment = ()
        self._load_abu(self._path(name, self.default_sets, silent=silent))
        self._apply_limits(**kwargs)
        self.is_sorted = True
        self.sort()

    def _path(self, name, sets, silent=False):
        """
        Find file from path or set name.

        TODO - should probably be named differently
        """
        self.setup_logger(silent=silent)
        if name in sets:
            self.name = name
            path_ = os.getenv("KEPLER_DATA")
            if not path_:
                path_ = path
                self.logger.warning(f"using default path {path_!s}")
            else:
                path_ = Path(path_)
            filename = path_ / sets[name]
        else:
            self.name = name
            filename = name
            if not os.path.isfile(filename):
                path_ = os.getenv("KEPLER_DATA")
                if not path_:
                    path_ = path
                    self.logger.warning(f"using default path {path_!s}")
                else:
                    path_ = Path(path_)
                filename = path_ / filename
        # HERE we should add error treatment
        if isinstance(filename, str):
            filename = Path(filename)
        if not filename.is_file():
            s = f"Data file for {filename} not found."
            self.logger.critical(s)
            raise Exception(f" [SolAbu] {s}")
        self.close_logger()
        return filename

    def _load_abu(self, filename):
        self._from_dat(filename, silent=self.silent)


def bbncoc(filename, write=False):
    """
    convert BBN abunaces from
    http://www2.iap.fr/users/pitrou/primat.htm
    to data file
    """

    with open(filename) as f:
        lines = f.readlines()
    ions = []
    abu = []
    for i, l in enumerate(lines):
        if l.count("1 n") > 0:
            break
    while True:
        if lines[i].count("Z=") > 0:
            break
        offset = 0
        if lines[i].startswith(";"):
            offset += 2
        xions = [
            "".join(lines[i][k : k + 10].split())
            for k in range(offset, 70 + offset, 10)
        ]
        i += 1
        xabu = list(lines[i][offset:].split())
        ions += [x for x in xions if len(x) > 0]
        abu += xabu
        i += 2
    abu = [float(a) for a in abu]
    ions = I(ions)
    abu = AbuSet(ions, abu)
    from ionmap import decay

    abu = decay(abu, stable=True)
    if write:
        with open(filename, "at") as f:
            for a in abu:
                f.write(f"{a[0].name():5s} {a[1]:12.5e}\n")
    return abu


def bbncyburt(filename=None):
    """
    BBN from Cyburt+2016
    https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.88.015004
    """
    Y = 0.24709
    D = 2.58e-5 * 2
    He3 = 10.039e-6 * 3
    Li7 = 4.68e-10 * 7
    Li6 = 10 ** (-13.89) * 6

    B = D + He3 + Li7 + Li6
    X = (1 - Y) / (1 + B)
    D *= X
    He3 *= X
    Li6 *= X
    Li7 *= X
    abu = [X, D, He3, Y, Li6, Li7]
    ions = I(["h1", "h2", "he3", "he4", "li6", "li7"])
    abu = AbuSet(ions, abu)
    if filename is not None:
        with open(filename, "at") as f:
            for a in abu:
                f.write(f"{a[0].name():5s} {a[1]:12.5e}\n")
    return abu


class BBNAbu(SolAbu):
    """
    load BBN set
    """

    default_sets = {
        "F02": "bbnf02.dat",
        "C16": "bbnc16.dat",
        "C19": "bbnc19.dat",
    }
    default = "C19"


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
            (jj,) = np.argwhere(scaled.iso == I.He4)
            # make sure we have same set of istopes
            bbn = (self.sun * 0) + self.bbn
            for j in np.argwhere(scaled.abu < self.sun.abu).flat:
                scaled.abu[jj] += scaled.abu[j]
                scaled.abu[j] = self.sun.abu[j] * np.exp(
                    (scale - 1) * (1 - bbn.abu[j] / self.sun.abu[j])
                )
                scaled.abu[jj] -= scaled.abu[j]
        scaled.normalize()

        return scaled.abu

    # use same defaults and definitions as in SolAbu
    def __init__(self, **kw):

        """
        Create abundance from set name.

        parameters:
          z:
            scale on solar abundance
          Z:
            alternatively, absolute metallicity (mass fraction)

        Use simple algorithm:
        X = X_sun * z + X_BBN * (1 - z)

        If Z is provided, overwrite z by Z/Zsun.

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

        simple = kw.pop("simple", False)
        if simple:
            Z = kw.pop("Z", None)
            z = kw.pop("z", None)
            if Z is not None:
                assert z is None
                z = -Z
            if z <= 0:
                Zsun = self.sun.XYZ()[-1]
                z = -z / Zsun
            self.z = z
            abu = self.sun * self.z + self.bbn * (1 - self.z)
            abu.normalize()
            AbuSet.__init__(self, abu)
            self.is_sorted = self.sun.is_sorted
        else:
            super().__init__(**kw)

        xyz = "X = {:8G}, Y = {:8G}, Z = {:8G}".format(*self.XYZ())
        self.comment = (
            f"Version {self.version:6s} - {time.asctime(time.gmtime())} UTC",
            f"Scaled solar abundances: {self.z} solar",
            f"Sun: {solar} - {self.sun.filename}",
            f"BBN: {zero} - {self.bbn.filename}",
            xyz,
        )
        self.logger.info(xyz)
        self.close_logger()


class ScaledSolarHelium(ScaledSolar):
    """
    Special abundance set created from scaled solar abundance set and
    overwrite He abundance.

    Version history:
    10000: created
    10001: add light isotope scaling
    10002: add fhelium and dhelium
    """

    version = "10002"

    def _abu_massfrac_raw(self, scale):
        """
        Raw scaled solar abundances.
        """
        # TODO - use abuset object for abu

        abu = super()._abu_massfrac_raw(scale)

        jH1, jHe4 = self.sun.index(("H1", "He4"))

        if self.dhelium is not None:
            self.helium = abu[jHe4] + self.dhelium
        elif self.fhelium is not None:
            self.helium = abu[jHe4] * self.fhelium

        if self.helium is None:
            return abu

        if self.scale_light is False:
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

            jH2, jHe3 = self.sun.index((I.H2, I.He3))
            # jLi6 = np.argwhere(self.sun.iso == I.Li6)
            # jLi7 = np.argwhere(self.sun.iso == I.Li7)
            # jBe9 = np.argwhere(self.sun.iso == I.Be9)
            # jB10 = np.argwhere(self.sun.iso == I.B10)
            # jB11 = np.argwhere(self.sun.iso == I.B11)

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
    def __init__(self, *args, **kwargs):

        (
            """
        Based on scaled solar, overwrite He4, rest in H1.

        parameters:
          z:
            scale on solar abundance
          Z:
            alternatively, absolute metallicity (mass fraction)

        (use one of the following three)
          dhelium:
            change helium by this mass fraction
          fhelium:
            change helium by this factor
          helium:
            desired He4 mass fraction specified explicitly

          scale_light:
            also scale other light istopes assumed to be
            destroyed to same fraction as H1

        Keeps H2 and He3 values by default.

        Base class documentation:
        -------------------------
        """
            + super().__doc__
        )

        self.helium = kwargs.pop("helium", None)
        self.dhelium = kwargs.pop("dhelium", None)
        self.fhelium = kwargs.pop("fhelium", None)
        if (
            np.count_nonzero(
                np.array([self.helium, self.dhelium, self.fhelium]) is not None
            )
            > 1
        ):
            raise AttributeError("Not more than one helium specification allowed.")

        self.scale_light = kwargs.pop("scale_light", False)
        super().__init__(*args, **kwargs)

        if self.helium is None:
            self.helium = self(I.He4)
            self.comment += ("Using scaled solar helium",)
        self.comment += (f"Helium set to mass fraction: {self.helium:g}",)


class _Asplund2009Data(AbuSet):
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
        comment = stuple(comment)

        xre = re.compile("[-+a-zA-Z0-9.]+")
        iso = np.array([], dtype=object)
        abu = np.array([], dtype=np.float64)

        filename = Path(filename).expanduser()

        with open(filename, "r") as f:
            self.logger_file_info(f)
            comment += (
                "",
                f'Generated from file "{filename}".',
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
                        xiso = I(
                            A=int(xdata[2]),
                            Z=int(xdata[1]),
                        )
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


def _Lodders2022():
    import abusets
    import bdat
    import ionmap
    from physconst import YR

    s = abusets.SolAbu("sollo20.dat")
    t = s.log_scaled_by(Y=-0.07, Z=-0.088)
    b = bdat.BDat("/home/alex/kepler/bdat_projjwal/bdat").decaydata
    d = ionmap.TimedDecay(ions=t, decaydata=b)
    tau = 4.5673e9 * YR
    u = d(t, tau, 1000)
    v = ionmap.decay(u, stable=True)
    f = v.Y("li") / v.Y("si") * 1e6
    ff = 0.339 / f
    ff1 = 1 - ff
    li = AbuSet(
        dict(
            h1=-v.Y("li") * ff1,
            he3=3 * v.Y("li6") * ff1,
            he4=4 * v.Y("li6") * ff1 + 8 * v.Y("li7") * ff1,
            li6=-v.li6 * ff1,
            li7=-v.li7 * ff1,
        ),
        allow_negative=True,
    )
    w = v + li
    w.write_dat(sys.stdout)


# seeems to be not used at present
# add option to AbuSet etc. to use proper masses
# use bdat
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

    # TODO - update to use different sets

    def __init__(self, data=default_data, silent=False):
        """
        TODO - implement data
        """
        self.setup_logger(silent=silent)
        path = os.getenv("KEPLER_DATA")
        if not path:
            path = Path("~").expanduser() / "kepler" / "local_data"
            self.logger.warning("using default path " + path)
        filename = path / "masses_audi_2003.dat"

        self.comment = ()
        self.iso = np.array([], dtype=object)
        self.mass = np.array([], dtype=np.float64)

        xre = re.compile("[-+a-zA-Z0-9.]+")
        with open(filename, "rt") as f:
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
                    self._append(I(xion), np.double(xabu))
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
        return np.double(I(ion).A)

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
        return ufunc_A(self.iso)

    @CachedAttribute
    def Z(self):
        return ufunc_Z(self.iso)

    @CachedAttribute
    def N(self):
        return ufunc_N(self.iso)

    @CachedAttribute
    def E(self):
        return ufunc_E(self.iso)

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
