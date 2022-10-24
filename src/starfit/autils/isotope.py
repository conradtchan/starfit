"""
Collection of routines for isotopes, most prominently the Ion class.

Usage:
------

use 'ion' factory function for many practical purposes, do not create
ions in other ways.

TODO - replace awful ion decpmposition by re matching


"""

# contains all the isotope routines

import collections
import re
import sys
from collections.abc import Iterable, Mapping

import numpy as np

# the following now is a singleton metaclass


class MetaSingletonHash(type):
    """
    Singleton metaclass based on hash

    First creates object to be able to test hash.

    If same hash is found, return old object and discard new one,
    otherwise return old one.

    Usage:
       class X(Y, metaclass = MetaSingletonHash)

    class X needs to provide a __hash__ function
    """

    def __call__(*args, **kwargs):
        cls = args[0]
        try:
            cache = cls._cache
        except:
            cache = dict()
            cls._cache = cache
        obj = type.__call__(*args, **kwargs)
        key = (cls.__name__, obj.__hash__())
        return cache.setdefault(key, obj)


# # maybe this should be an enumeration???
# TODO - this class draft needs more work for tu
# class _elements(object):
#     _data = (
#     'nt',  'h',   'he',  'li',  'be',  'b',   'c',   'n',
#     'o',   'f',   'ne',  'na',  'mg',  'al',  'si',  'p',
#     's',   'cl',  'ar',  'k',   'ca',  'sc',  'ti',  'v',
#     'cr',  'mn',  'fe',  'co',  'ni',  'cu',  'zn',  'ga',
#     'ge',  'as',  'se',  'br',  'kr',  'rb',  'sr',  'y',
#     'zr',  'nb',  'mo',  'tc',  'ru',  'rh',  'pd',  'ag',
#     'cd',  'in',  'sn',  'sb',  'te',  'i',   'xe',  'cs',
#     'ba',  'la',  'ce',  'pr',  'nd',  'pm',  'sm',  'eu',
#     'gd',  'tb',  'dy',  'ho',  'er',  'tm',  'yb',  'lu',
#     'hf',  'ta',   'w',  're',  'os',  'ir',  'pt',  'au',
#     'hg',  'tl',  'pb',  'bi',  'po',  'at',  'rn',  'fr',
#     'ra',  'ac',  'th',  'pa',  'u',   'np',  'pu',  'am',
#     'cm',  'bk',  'cf',  'es',  'fm',  'md',  'no',  'lr',
#     'rf',  'db',  'sg',  'bh',  'hs',  'mt',  'ds',  'rg',
#     'cn',  'ut',  'fl',  'up',  'lv',  'us',  'uo')

#     def __getitem__(self, index):
#         if isinstance(index, str):
#             try:
#                 index = int(index)
#             except:
#                 pass
#         try:
#             return self._data[index]
#         except:
#             pass
#         try:
#             return self.index(index)
#         except:
#             raise IndexError()
#     def __len__(self):
#         return len(self._data)
#     def index(self, item):
#         return self._data.index(item.lower())
#     def __call__(self, value):
#         return self.__getitem__(value)
#     def __getattr__(self, attr):
#         try:
#             return self.index(attr)
#         except:
#             raise AttributeError()
# elements = _elements()

# should be Uus Uuo but I will not use 3 letters
# as this breaks too many things

# class _Elements(_elements):
#     _data = [e.capitalize() for e in _elements._data]
#     _data[0] = 'nt'
#     _data = tuple(_data)
#     def index(self, value):
#         try:
#             return self._data[value]
#         except ValueError:
#             return self._data[value.capitalize()]
# Elements = _Elements()

elements = (
    "nt",
    "h",
    "he",
    "li",
    "be",
    "b",
    "c",
    "n",
    "o",
    "f",
    "ne",
    "na",
    "mg",
    "al",
    "si",
    "p",
    "s",
    "cl",
    "ar",
    "k",
    "ca",
    "sc",
    "ti",
    "v",
    "cr",
    "mn",
    "fe",
    "co",
    "ni",
    "cu",
    "zn",
    "ga",
    "ge",
    "as",
    "se",
    "br",
    "kr",
    "rb",
    "sr",
    "y",
    "zr",
    "nb",
    "mo",
    "tc",
    "ru",
    "rh",
    "pd",
    "ag",
    "cd",
    "in",
    "sn",
    "sb",
    "te",
    "i",
    "xe",
    "cs",
    "ba",
    "la",
    "ce",
    "pr",
    "nd",
    "pm",
    "sm",
    "eu",
    "gd",
    "tb",
    "dy",
    "ho",
    "er",
    "tm",
    "yb",
    "lu",
    "hf",
    "ta",
    "w",
    "re",
    "os",
    "ir",
    "pt",
    "au",
    "hg",
    "tl",
    "pb",
    "bi",
    "po",
    "at",
    "rn",
    "fr",
    "ra",
    "ac",
    "th",
    "pa",
    "u",
    "np",
    "pu",
    "am",
    "cm",
    "bk",
    "cf",
    "es",
    "fm",
    "md",
    "no",
    "lr",
    "rf",
    "db",
    "sg",
    "bh",
    "hs",
    "mt",
    "ds",
    "rg",
    "cn",
    "nh",
    "fl",
    "mc",
    "lv",
    "ts",
    "og",
    "ue",
)

Elements = [element.capitalize() for element in elements]
Elements[0] = "nt"
Elements = tuple(Elements)

z2Name = (
    "Neutron",
    "Hydrogen",
    "Helium",
    "Lithium",
    "Beryllium",
    "Boron",
    "Carbon",
    "Nitrogen",
    "Oxygen",
    "Fluorine",
    "Neon",
    "Sodium",
    "Magnesium",
    "Aluminium",
    "Silicon",
    "Phosphorus",
    "Sulfur",
    "Chlorine",
    "Argon",
    "Potassium",
    "Calcium",
    "Scandium",
    "Titanium",
    "Vanadium",
    "Chromium",
    "Manganese",
    "Iron",
    "Cobalt",
    "Nickel",
    "Copper",
    "Zinc",
    "Gallium",
    "Germanium",
    "Arsenic",
    "Selenium",
    "Bromine",
    "Krypton",
    "Rubidium",
    "Strontium",
    "Yttrium",
    "Zirconium",
    "Niobium",
    "Molybdenum",
    "Technetium",
    "Ruthenium",
    "Rhodium",
    "Palladium",
    "Silver",
    "Cadmium",
    "Indium",
    "Tin",
    "Antimony",
    "Tellurium",
    "Iodine",
    "Xenon",
    "Caesium",
    "Barium",
    "Lanthanum",
    "Cerium",
    "Praseodymium",
    "Neodymium",
    "Promethium",
    "Samarium",
    "Europium",
    "Gadolinium",
    "Terbium",
    "Dysprosium",
    "Holmium",
    "Erbium",
    "Thulium",
    "Ytterbium",
    "Lutetium",
    "Hafnium",
    "Tantalum",
    "Tungsten",
    "Rhenium",
    "Osmium",
    "Iridium",
    "Platinum",
    "Gold",
    "Mercury",
    "Thallium",
    "Lead",
    "Bismuth",
    "Polonium",
    "Astatine",
    "Radon",
    "Francium",
    "Radium",
    "Actinium",
    "Thorium",
    "Protactinium",
    "Uranium",
    "Neptunium",
    "Plutonium",
    "Americium",
    "Curium",
    "Berkelium",
    "Californium",
    "Einsteinium",
    "Fermium",
    "Mendelevium",
    "Nobelium",
    "Lawrencium",
    "Rutherfordium",
    "Dubnium",
    "Seaborgium",
    "Bohrium",
    "Hassium",
    "Meitnerium",
    "Darmstadtium",
    "Roentgenium",
    "Copernicium",
    "Nihonium",
    "Flerovium",
    "Moscovium",
    "Livermorium",
    "Tennessine",
    "Oganesson",
    "(Ununennium)",
)
z2name = tuple((name.lower() for name in z2Name))


# el2z provides easy conversion from element symbol to charge number
el2z = {el: i for i, el in enumerate(elements)}
el2z.update({el: i for i, el in enumerate(Elements[1:], start=1)})
el2z["d"] = 1
el2z["t"] = 1
el2z["D"] = 1
el2z["T"] = 1

elements_ = elements + ("t", "d", "T", "D")

iso_special_mass = dict(
    d=2,
    D=2,
    T=3,
    t=3,
)
iso_special_element = ("D", "T")
iso_special_particle = ("d", "t")

iso_strings_ground = ("g", "G")
iso_strings_excite = ("m", "M")

iso_string_ground = "g"
iso_string_excite = "m"

iso_string_eany = "*"

iso_strings_ = ("g", "m")
iso_strings = ("g", "m", "G", "M")

iso_strings_all = iso_strings + (iso_string_eany,)


def el2name(sym):
    return z2name[el2z[sym]]


class Ion(object, metaclass=MetaSingletonHash):
    """
    This is a generic class for 'ions'.
    It is designed to handle elements, isotioes, isomers, and isobars.
    Isomers have an energy level appended by 'e' or 'E';
        if the level is '0' then it should not be printed.
    It provides also an "index" field (idx) that can be used for sorting.
    __add__ and __sub__ should help with networks.
    R([[I],O]) allows to formulate reactions
        will return ground state nuclei
    Use (solely)
        A = ... for 'isobars'
        Z = ... for 'elements'
        N = ... for 'isotones'
    To specify isotopes you need to supply at least 2 of A, Z, N;
        if you specify three, they must fulfil A = N + Z
    To specify an isomer, you have provide additionally E = ...
        Use E = 0 to specify an ismer in the gs.

    g    = ( 0, 0, 0, 1)
    e-   = ( 7,-1, 0, 1)
    e+   = ( 7,+1, 0, 2)
    void = (-1, 0, 0, 0)

    We always put E=1, i.e., m1, as 'm'
    Default use 'mE' to indicate exited state E.
    We distinguish isotopes and isomers - Al26g is not Al26.
    Use EANY to specify generalised excited state, signified by '*'.
    """

    # in the future we want to use these
    # should flags be first or last?
    # maybe last?
    AMAX = 512
    ZMAX = 128  # -7 ... 120
    ZLIM = 120
    EMAX = 256
    FMAX = 128  # flags - 7 bit
    # leave the sign bit for int64

    FMUL = 1
    EMUL = FMUL * FMAX
    AMUL = EMUL * EMAX
    ZMUL = AMUL * AMAX
    BMUL = ZMUL * ZMAX
    # what lies Beyond ... (extra information for subclasses)
    # if these are used, the index no longer fits into int64

    # bit 0-2
    F_BOSON = 0
    F_ELEMENT = 1
    F_ISOTOPE = 2
    F_ISOMER = 3
    F_ISOBAR = 4
    F_ISOTONE = 5
    F_IONIZE = 6  # ionisation stage (element or isotope
    F_LEPTON = 7
    F_GROUP_MASK = 7

    # bits 3-6 remain for other implementation
    F_OTHER = 8
    F_OTHER_MASK = F_OTHER * (FMAX // F_OTHER - 1)

    _Excited = "E"
    _excited = "e"
    ISOMER = "m"
    GROUND = "g"
    EXCITE = (_Excited, _excited, ISOMER)
    EANY = EMAX - 1
    EXCITED_ANY = "*"
    VOID = (0, -1, 0, 0, 0)
    VOID_IDX = -1 * FMUL
    VOID_STRING = "-"
    GAMMA = (0, F_BOSON, 0, 0, +1)
    VAC = (0, 0, 0, 0, 0)
    VAC_STRING = "."
    VAC_IDX = 0

    E_LEPTON_FLAVOR_MUL = 2
    # lepton number is 2 * (E % 2) - 1, stored a 0, 1 for -1, 1
    # lepton flavors are returned 1-3, stored as 0, 2, 4

    E_LEPTON = 1
    E_LEPTON_ANTI = 0
    E_LEPTON_E = 1
    E_LEPTON_M = 2
    E_LEPTON_T = 3
    E_LEPTON_X = 4
    E_LEPTON_ANY = EMAX // 2

    # E = E_LEPTON_ANY|(E_LEPTON|E_LEPTON_ANTI)+(E_LEPTON_E|E_LEPTON_M|E_LEPTON_T)*E_LEPTON_FLAVOR_MUL

    _SPECIAL = {
        0: "nt",
        1 * EMUL: "g",
        (F_LEPTON * FMUL) + (+0 * EMUL) + (+1 * ZMUL): "e+",
        (F_LEPTON * FMUL) + (+1 * EMUL) + ((ZMAX - 1) * ZMUL): "e-",
        (F_LEPTON * FMUL) + (+2 * EMUL) + (+1 * ZMUL): "m+",
        (F_LEPTON * FMUL) + (+3 * EMUL) + ((ZMAX - 1) * ZMUL): "m-",
        (F_LEPTON * FMUL) + (+4 * EMUL) + (+1 * ZMUL): "t+",
        (F_LEPTON * FMUL) + (+5 * EMUL) + ((ZMAX - 1) * ZMUL): "t-",
        (F_LEPTON * FMUL) + (+0 * EMUL): "nbe",
        (F_LEPTON * FMUL) + (+1 * EMUL): "nue",
        (F_LEPTON * FMUL) + (+2 * EMUL): "nbm",
        (F_LEPTON * FMUL) + (+3 * EMUL): "num",
        (F_LEPTON * FMUL) + (+4 * EMUL): "nbt",
        (F_LEPTON * FMUL) + (+5 * EMUL): "nut",
        (F_LEPTON * FMUL)
        + ((E_LEPTON_X - 1) * 2 * EMUL): "nbx",  # generic anti-neutrino
        (F_LEPTON * FMUL)
        + (((E_LEPTON_X - 1) * 2 + 1) * EMUL): "nux",  # generic lepton
        (F_LEPTON * FMUL)
        + ((E_LEPTON_X - 1) * 2 * EMUL)
        + (+1 * ZMUL): "l+",  # generic anti-lepton
        (F_LEPTON * FMUL)
        + (((E_LEPTON_X - 1) * 2 + 1) * EMUL)
        + ((ZMAX - 1) * ZMUL): "l-",  # generic lepton
        (F_LEPTON * FMUL) + ((E_LEPTON_ANY - 1) * 2 * EMUL): "nu",  # generic neutrino
        (F_LEPTON * FMUL)
        + (((E_LEPTON_ANY - 1) * 2 + 1) * EMUL): "l",  # generic lepton
        VOID_IDX: VOID_STRING,
        VAC_IDX: VAC_STRING,
    }
    SPECIAL = {v: k for k, v in _SPECIAL.items()}

    def __init__(
        self,
        name=None,
        Z=None,
        A=None,
        E=None,
        F=None,
        B=None,
        idx=None,
        N=None,
        element=False,
        isomer=False,
        isobar=False,
        isotone=False,
        lepton=False,
        ionize=False,
        void=True,
        factory=False,
        check=None,
        _load=None,
    ):
        """
        Initialize Isotope.

        isomer == True:
            interpret istopes as gs isomers
        """

        if _load is not None:
            try:
                # new version
                self.B, self.F, self.Z, self.A, self.E, self.N, self.idx = _load
            except:
                self.B, self.F, self.Z, self.A, self.E = _load
                self.N = self._get_N()
                self.idx = self.bfzae2idx(
                    B=self.B, F=self.F, Z=self.Z, A=self.A, E=self.E
                )
            return

        if self.__class__ == Ion and factory is False:
            raise NotImplementedError("You should not call Ion constructor directly")

        if check is None:
            check = isinstance(name, str)

        # experimental
        if A == -1:
            A = None
        if N == -1:
            N = None
        if E == -1:
            E = None
        if not lepton:
            if Z == -1:
                Z = None
        else:
            F = self.F_LEPTON if F is None else F

        self.F = self.F_ISOMER if F is None else F
        self.Z = 0 if Z is None else Z
        self.A = 0 if A is None else A
        self.E = 0 if E is None else E
        self.B = 0 if B is None else B

        if N is not None:
            if Z is not None:
                self.A = Z + N
                assert (A == self.A) or (A is None), "Conflict."
            elif A is not None:
                self.Z = A - N
                self.A = A
            else:
                self.F = self.F_ISOTONE
                self.A = N
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
                assert self.A >= 0, "Not sure what this would be."
        if Z is not None:
            if (N is None) and (A is None) and (E is None):
                self.F = self.F_ELEMENT
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
            elif B is None and not lepton and not self.F == self.F_IONIZE:
                assert self.A != 0, "Not sure what this would be."
        if A is not None:
            if (Z is None) and (N is None):
                self.F = self.F_ISOBAR
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
                assert self.A > 0, "Not sure what this would be."
        if (E is None) or (E == -1):
            if self.F & self.F_GROUP_MASK == self.F_ISOMER:
                if not isomer:
                    self.F = self.F_ISOTOPE
                self.E = 0

        s = name
        if self._is_ion(s):
            self.B, self.F, self.Z, self.A, self.E = s.tuple
            # conversion to Ion class
            # or could add keyword
            if not issubclass(type(self), type(s)):
                self.B = 0
                self.F = self.F & self.F_GROUP_MASK
            s = None
            check = False

        if idx is None:
            if isinstance(s, str):
                try:
                    s = complex(s)
                except ValueError:
                    pass
            try:
                assert np.isscalar(s)
                assert not np.iscomplex(s)
                i = int(np.real(s))
            except (ValueError, TypeError, AssertionError):
                pass
            else:
                if isobar:
                    self.A = i
                    assert self.Z == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ISOBAR
                    assert (self.F == F) or (F is None), "Conflict."
                elif isotone:
                    self.A = i
                    assert self.Z == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ISOTONE
                    assert (self.F == F) or (F is None), "Conflict."
                elif element:
                    self.Z = i
                    assert self.A == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ELEMENT
                    assert (self.F == F) or (F is None), "Conflict."
                elif lepton:
                    self.Z = i
                    if self.Z == 0 and E is None:
                        self.E = 1
                    else:
                        self.E = E if E is not None else 0
                    assert self.A == 0
                    assert self.B == 0
                    self.F = self.F_LEPTON
                    assert (self.F == F) or (F is None), "Conflict."
                elif 0 < i < len(Elements):
                    self.Z = i
                    self.A = Z if Z is not None else 0
                    self.E = A if A is not None else 0
                    self.F = self.F_ISOMER
                    if A is None:
                        self.F = self.F_ISOTOPE
                    if Z is None:
                        self.F = self.F_ELEMENT
                    assert (self.F == F) or (F is None), "Conflict."
                else:
                    idx = s
                s = None

        if isinstance(s, str):
            try:
                x = tuple(eval(re.sub(r"[- \.;:/\|]", ",", s)))
                s = tuple(int(i) for i in x)
            except:
                pass

        if isinstance(s, tuple):
            if len(s) == 0:
                self.F = 0
            if len(s) == 1:
                self.Z = s[0]
                self.F = self.F_ELEMENT
            if len(s) == 2:
                self.Z, self.A = s
                self.F = self.F_ISOTOPE
            if len(s) == 3:
                self.Z, self.A, self.E = s
                self.F = self.F_ISOMER
            if len(s) == 4:
                self.F, self.Z, self.A, self.E = s
            if len(s) == 5:
                self.B, self.F, self.Z, self.A, self.E = s
            s = None

        if np.iscomplex(s):
            assert Z is None
            assert A is None
            assert N is None
            self.Z = np.real(s)
            x = np.imag(s)
            s = None
            if x > 0:
                self.A = x
            else:
                self.A = self.Z - x
            if E is not None:
                self.E = E
                self.F = self.F_ISOMER
            else:
                self.F = self.F_ISOTOPE

        if isinstance(s, Ion):
            self.B, self.F, self.Z, self.A, self.E = s.tuple
            check = False
        elif s is not None:
            self.B, self.F, self.Z, self.A, self.E = self.ion2bfzae(
                s, element=element, isomer=isomer
            )
        if idx is not None:
            self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(idx)
        self.B = int(self.B)
        self.F = int(self.F)
        self.Z = int(self.Z)
        self.A = int(self.A)
        self.E = int(self.E)
        self.idx = self.bfzae2idx(B=self.B, F=self.F, Z=self.Z, A=self.A, E=self.E)
        # some checks?
        if self.F == self.F_IONIZE:
            if self.E > self.Z or ((self.A > 0) and (self.A < self.Z)):
                self.idx = self.VOID_IDX
                self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(self.idx)
            if self.A > 0 and check:
                if not self._valid_isotope(self.A, self.Z):
                    self.idx = self.VOID_IDX
                    self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(self.idx)
        if self.F in (self.F_ISOTOPE, self.F_ISOMER) and check:
            if not self._valid_isotope(self.A, self.Z):
                self.idx = self.VOID_IDX
                self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(self.idx)
        if self.idx == self.VOID_IDX and not void:
            raise AttributeError("Could not create valid isotope from provided input.")
        self.N = self._get_N()
        if isinstance(self, Ion):
            super().__init__()

    _limits = {
        8: [12, 34],
        9: [14, 38],
        10: [16, 41],
        11: [18, 44],
        12: [19, 47],
        13: [21, 51],
        14: [22, 54],
        15: [23, 57],
        16: [24, 60],
        17: [25, 63],
        18: [27, 67],
    }

    @classmethod
    def _valid_isotope(cls, A, Z):
        if A < Z:
            return False
        if Z == 0 and A > 1:
            return False
        if Z in cls._limits:
            return A >= cls._limits[Z][0] and A <= cls._limits[Z][1]
        if A > 3.1 * Z + 8 + 0.01 * Z**2:
            return False  # B21, Na33
        if A < 1.6 * Z - 3:
            return False
        return True

    @staticmethod
    def _is_ion(ix):
        return isinstance(ix, Ion)

    @classmethod
    def min_iso(cls, Z):
        if issubclass(type(Z), cls):
            Z = Z.Z
        try:
            Z = int(Z)
        except:
            Z = el2z[Z]
        return cls(Z=Z, N=0)

    @classmethod
    def max_iso(cls, Z):
        if issubclass(type(Z), cls):
            Z = Z.Z
        try:
            Z = int(Z)
        except:
            Z = el2z[Z]
        return cls(Z=Z, A=cls.AMAX - 1)

    @classmethod
    def element_slice(cls, Z):
        return slice(cls.min_iso(Z), cls.max_iso(Z))

    @classmethod
    def bfzae2idx(cls, F=F_ISOTOPE, Z=0, A=0, E=0, B=0):
        """
        Compute index from F, Z, A, E.
        """
        assert (
            ((F == -1) and (A == 0) and (Z == 0) and (E == 0) and (B == 0))
            or ((F >= 0) and (A >= 0) and (E >= 0) and (B >= 0))
            or (F >= cls.F_OTHER)
        ), "Wrong numbers."
        if Z < 0:
            Z += cls.ZMAX
        return F * cls.FMUL + E * cls.EMUL + Z * cls.ZMUL + A * cls.AMUL + B * cls.BMUL

    @classmethod
    def idx2bfzae(cls, idx):
        """
        Return things that are defined within basic Ion class.
        """
        aidx = abs(idx)
        f = (aidx // cls.FMUL) % cls.FMAX
        z = (aidx // cls.ZMUL) % cls.ZMAX
        a = (aidx // cls.AMUL) % cls.AMAX
        e = (aidx // cls.EMUL) % cls.EMAX
        b = aidx // cls.BMUL

        if z > cls.ZLIM:
            z -= cls.ZMAX
        elif idx == cls.VOID_IDX:
            b, f, z, a, e = cls.VOID
        elif idx == cls.VAC_IDX:
            b, f, z, a, e = cls.VAC
        return b, f, z, a, e

    @classmethod
    def ion2idx(cls, s="", element=False, isomeer=False):
        return cls.bfzae2idx(cls.ion2bfzae(s, element=element, isomer=isomer))

    def _update_idx(self):
        """
        Update idx assuming B,F,A,Z,E are all current
        """
        self.idx = (
            self.F * self.FMUL
            + self.E * self.EMUL
            + self.Z * self.ZMUL
            + self.A * self.AMUL
            + self.B * self.BMUL
        )

    def _get_N(self):
        N = self.A - self.Z
        if self.F in (
            self.F_BOSON,
            self.F_ELEMENT,
            self.F_ISOBAR,
            self.F_IONIZE,
            self.F_LEPTON,
        ):
            N = 0
        return N

    @property
    def tuple(self):
        return (self.B, self.F, self.Z, self.A, self.E)

    @property
    def dict(self):
        return {"F": self.F, "Z": self.Z, "A": self.A, "E": self.E, "B": self.B}

    @property
    def index(self):
        return self.idx

    @property
    def index_ZAE(self):
        return self.E * self.EMUL + self.Z * self.ZMUL + self.A * self.AMUL

    @property
    def ZAE(self):
        return (
            self.Z,
            self.A,
            self.E,
        )

    # A,Z, ... should become properties that are set on init

    @property
    def Ye(self):
        if self.Z < 0 or self.A == 0:
            return np.nan
        return self.Z / self.A

    ye = Ye

    @property
    def mu(self):
        if self.Z < 0:
            return np.nan
        return self.A / (self.Z + 1)

    @property
    def mue(self):
        if self.Z <= 0:
            return np.nan
        return self.A / self.Z

    @property
    def eta(self):
        return 1 - 2 * self.Ye

    def Name(self, align=1, width=0, upcase=True):
        """
        Return capitalized ion name; same as name otherwise.
        """
        return self.name(align=align, width=width, upcase=upcase)

    def name(self, align=1, width=0, upcase=False):
        """
        Return ion name.

        align [+1]:
            -1: left
             0: center
            +1: right (default)
        width [0]:
            size of field
        upcase [False]:
            capitalize first letter
        """
        s = self._name(upcase=upcase)
        if width > 0:
            if align == -1:
                s = s.ljust(width)
            elif align == 1:
                s = s.rjust(width)
            else:
                s = s.center(width)
        return s

    def element_symbol(self, upcase=True):
        s = ""
        if self.Z >= 0:
            if upcase:
                s = Elements[self.Z]
            else:
                s = elements[self.Z]
        return s

    def _math(func, math_=False):
        def f(*args, **kwargs):
            math = kwargs.pop("math", math_)
            s = func(*args, **kwargs)
            if math:
                return s
            return "$" + s + "$"

        return f

    def _LaTeX(self, **kwargs):

        if self.is_photon:
            s = r"\gamma"
            if self.E > 1:
                s = rf"{{{s}}}^{{{self.E}}}"
            return s

        if self.is_lepton:
            f = (r"\mathrm{e}", r"\mu", r"\tau", r"\mathrm{l}")[self.E // 2]
            if self.Z == 0:
                if self.E % 2 == 1:
                    l = r"\nu"
                else:
                    l = r"\bar{\nu}"
                s = rf"{{{l}}}_{{{f}}}"
            else:
                if self.Z == 1:
                    c = "+"
                else:
                    c = "-"
                s = rf"{{{f}}}^{{{c}}}"
            return s

        el = Elements[self.Z]
        if self.is_element or (self.is_ionize and self.A == 0):
            if self.Z == 0:
                el = "n"
            s = rf"\mathrm{{{el:s}}}"
            if self.is_element:
                return s

        if self == NEUTRON:
            return r"\mathrm{n}"

        a = str(self.A)
        if self.is_isotope or (self.is_ionize and self.A > 0):
            if el.startswith("A"):
                el = r"\!" + el
            elif el.startswith("C"):
                el = r"\!" + el
            s = rf"^{{{a:s}}}\mathrm{{{el:s}}}"
            if self.is_isotope:
                return s

        if self.is_ionize:
            i = self.ionize_name(self.E)
            return s + rf"\,\mathsc{{{i}}}"

        if self.is_isomer:
            kw = {k: kwargs[k] for k in ("m1", "g") if k in kwargs}
            e = rf"\mathrm{{{self.isomer_name(self.E, **kw)}}}"
            return rf"^{{{a:s}}}\mathrm{{{el:s}}}^{{{e:s}}}"

        if self.is_isobar:
            return rf"A\!\!=\!\!{self.A}"

        return "-"

    LaTeX = _math(_LaTeX, math_=False)
    latex = _math(_LaTeX, math_=True)

    @property
    def NuGrid(self):
        el = Elements[self.Z]
        if self.is_element:
            return el
        A = str(self.A)
        if self.is_isotope:
            return f"{el}-{A}"
        if self.is_isomer:
            return f"{el}-{A}{self.isomer_name(self.E, m1=True)}"
        raise Exception("Cannot convert ion to string.")

    @property
    def nugrid5(self):
        el = elements[self.Z]
        if self.is_isotope:
            s = f"{el:2s}{self.A:3d}"
            if s == "nt  1":
                s = "NEUT"
            elif s == "h   1":
                s = "PROT"
            return s
        raise Exception("Cannot convert ion to string.")

    @property
    def nugrid5z(self):
        if self.is_isotope:
            return f"{self.Z:3d} {self.nugrid5:5s}"
        raise Exception("Cannot convert ion to string.")

    @property
    def Kepler(self):
        el = elements[self.Z]
        if self.is_element:
            return el
        A = str(self.A)
        if self.is_isotope:
            return f"{el}{A}"
        if self.is_isomer:
            return f"{el}{A}{self.isomer_name(self.E, m1=True)}"
        raise Exception("Cannot convert ion to string.")

    @classmethod
    def isomer_name(cls, E, m1=False, g=True):
        """
        For isomers use 'g', 'm', 'm2', m3', ...
        """
        if E == 0:
            if g is False:
                return ""
            return cls.GROUND
        elif E == 1 and m1 is False:
            return f"{cls.ISOMER:s}"
        elif E == cls.EANY:
            return f"{cls.EXCITED_ANY:s}"
        else:
            return f"{cls.ISOMER:s}{E:d}"

    @property
    def mpl(self):
        if self.is_element:
            el = Elements[self.Z]
            return rf"$\mathsf{{{el}}}$"
        if self.is_isobar:
            return rf"$\mathsf{{{self.A:d}}}$"
        if self.is_isotone:
            return rf"$\mathsf{{{self.N:d}}}$"
        if self.is_isotope:
            if self == NEUTRON:
                return r"$\mathsf{n}$"
            el = Elements[self.Z]
            a = str(self.A)
            return rf"$\mathsf{{^{{{a}}}{el}}}$"
        if self.is_isomer:
            el = Elements[self.Z]
            a = str(self.A)
            m = self.isomer_name(self.E)
            return rf"$\mathsf{{^{{{a}}}{el}^{{{m}\!}}}}$"
        if self.is_ionize:
            el = Elements[self.Z]
            i = self.ionize_name(self.E, lower=False)
            if self.A > 0:
                a = str(self.A)
                return rf"$\mathsf{{^{{{a}}}{el}\,{i}}}$"
            else:
                return rf"$\mathsf{{{el}\,{i}}}$"
        if self.is_lepton:
            f = (r"e", r"\mu", r"\tau", r"l")[self.E // 2]
            if self.Z == 0:
                if self.E % 2 == 1:
                    l = r"\nu"
                else:
                    l = r"\bar{\nu}"
                s = rf"{{{l}}}_{{{f}}}"
            else:
                if self.Z == 1:
                    c = "+"
                else:
                    c = "-"
                s = rf"{{{f}}}^{{{c}}}"
            return rf"$\mathsf{{{s}}}$"
        if self.is_photon:
            s = r"\gamma"
            if self.E > 1:
                s = rf"{{{s}}}^{{{self.E}}}"
            return rf"$\mathsf{{{s}}}$"
        return "#"

    @property
    def mpl_(self):
        return self.mpl.strip("$")

    def LaTeX_Table(
        self, math=False, nozero=False, protons=False, manual=False, width=11, align=0
    ):
        r"""
        Add the following definitions:

        \newcommand{\Z}[3][]{$^{#1}\mathsf{#2}\,\textsc{#3}$}
        \newcommand{\I}[3][]{$^{#2}\mathsf{#3}^{#1}$}
        """

        el = Element[self.Z]
        if manual:
            if self.is_isotope:
                a = str(self.A)
                s = f"$^{{{a}}}$&$\\mathsf{{{el}}}$"
            elif self.is_isomer:
                a = str(self.A)
                m = self.isomer_name(self.E)
                s = f"$^{{{a}}}$&$\\mathsf{{{el}}}$&$^{{{m}}}"
            elif self.is_ionize:
                i = self.ionize_name(self.E)
                if self.A > 0:
                    a = str(self.A)
                    s = f"$^{{{a}}}$&$\\mathsf{{{el}}}$&\\textsc{{{i}}}"
                else:
                    s = f"$\\mathsf{{{el}}}$&\\textsc{{{i}}}"
            elif self.is_element:
                s = el
        else:
            if self.is_isotope:
                a = str(self.A)
                s = f"\\I{{{a}}}{{{el}}}"
            elif self.is_isomer:
                a = str(self.A)
                m = self.isomer_name(self.E)
                s = f"\\I[{m}]{{{a}}}{{{el}}}"
            elif self.is_ionize:
                i = self.ionize_name(self.E)
                if self.A > 0:
                    a = str(self.A)
                    s = f"\\Z[{a}]{{{el}}}{{{i}}}"
                else:
                    s = f"\\Z{{{el}}}{{{i}}}"
            elif self.is_element:
                s = el

        if nozero and self.idx == self.VOID_IDX:
            s = "--"
        elif nozero and self.idx == self.VAC_IDX:
            s = self.VAC_STRING
        if align == -1:
            s = s.ljust(width)
        elif align == 1:
            s = s.rjust(width)
        else:
            s = s.center(width)
        if math and not manual:
            s = "$" + s + "$"
        return s

    _roman_numerals = (
        "I",
        "II",
        "III",
        "IV",
        "V",
        "VI",
        "VII",
        "VIII",
        "IX",
        "X",
        "XI",
        "XII",
        "XIII",
        "XIV",
        "XV",
        "XVI",
        "XVII",
        "XVIII",
        "XIX",
        "XX",
        "XXI",
        "XXII",
        "XXIII",
        "XXIV",
        "XXV",
        "XXVI",
        "XXVII",
        "XXVIII",
        "XXIX",
        "XXX",
        "XXXI",
        "XXXII",
        "XXXIII",
        "XXXIV",
        "XXXV",
        "XXXVI",
        "XXXVII",
        "XXXVIII",
        "XXXIX",
        "XV",
        "XLI",
        "XLII",
        "XLIII",
        "XLIV",
        "XLV",
        "XLVI",
        "XLVII",
        "XLVIII",
        "XLIX",
        "L",
        "LI",
        "LII",
        "LIII",
        "LIV",
        "LV",
        "LVI",
        "LVII",
        "LVIII",
        "LIX",
        "LX",
        "LXI",
        "LXII",
        "LXIII",
        "LXIV",
        "LXV",
        "LXVI",
        "LXVII",
        "LXVIII",
        "LXIX",
        "LXX",
        "LXXI",
        "LXXII",
        "LXXIII",
        "LXXIV",
        "LXXV",
        "LXXVI",
        "LXXVII",
        "LXXVIII",
        "LXXIX",
        "LXXX",
        "XCI",
        "XCII",
        "XCIII",
        "XCIV",
        "XCV",
        "XCVI",
        "XCVII",
        "XCVIII",
        "XCIX",
        "C",
        "CI",
        "CII",
        "CIII",
        "CIV",
        "CV",
        "CVI",
        "CVII",
        "CVIII",
        "CIX",
        "CX",
        "CXI",
        "CXII",
        "CXIII",
        "CXIV",
        "CXV",
        "CXVI",
        "CXVII",
        "CXVIII",
        "CXIX",
        "CXX",
    )

    def ionize_name(self, number=0, roman=True, lower=False):
        if roman:
            s = self._roman_numerals[number]
            if lower:
                s = s.lower()
            return s
        return f"({number+1:d})"

    def _name(self, upcase=True):
        """
        Convert normal isotopes and isomers to string.

        For isomers use 'g', 'm', 'm2', 'm3', ...

        For ionize use 'I', 'II', 'III', 'IV', ...
        """
        s = ""
        if self.Z >= 0:
            if upcase:
                s = Elements[self.Z]
            else:
                s = elements[self.Z]
        if self.F & self.F_GROUP_MASK == self.F_ISOBAR:
            s = "A:"
        if self.F & self.F_GROUP_MASK == self.F_ISOTONE:
            s = "N:"
        if self.A != 0 or (self.F & self.F_GROUP_MASK == self.F_ISOTONE):
            s += f"{self.A:d}"
        if (
            self.A == 1
            and self.Z == 0
            and (self.F & self.F_GROUP_MASK == self.F_ISOTOPE)
        ):
            s = "n"
        if self.F & self.F_GROUP_MASK == self.F_ISOMER:
            if self.A == 0 and self.Z == 0:
                if self.E == 1:
                    s = "g"
                else:
                    s = f"g{self.E:d}"
            else:
                s += self.isomer_name(self.E)
        if self.F & self.F_GROUP_MASK == self.F_BOSON:
            if self.A == 0 and self.Z == 0:
                if self.E == 1:
                    s = "g"
                else:
                    s = f"g{self.E:d}"
            else:
                raise NotImplementedError()
        if self.F & self.F_GROUP_MASK == self.F_IONIZE:
            if self.A == 0 and self.Z == 0:
                if self.E == 1:
                    s = "g"
                else:
                    s = f"g{self.E:d}"
            else:
                s += self.ionize_name(self.E)
        s = self._SPECIAL.get(self.idx, s)
        return s

    @property
    def lepton_number(self):
        """
        return lepton number
        """
        if self.is_lepton:
            if self.E // 2 == self.E_LEPTON_ANY - 1:
                return 0
            return 2 * (self.E % 2) - 1
        return 0

    @property
    def lepton_flavor(self):
        """return lepton flavor 1-3, 0 for non-leptons"""
        if self.is_lepton:
            if self.E // 2 == self.E_LEPTON_ANY - 1:
                return 0
            return self.E // 2 + 1
        return 0

    @property
    def lepton_flavor_vec(self):
        """return lepton flavor vector [l_e, l_nu, l_tau]"""
        v = [0, 0, 0]
        if self.is_lepton:
            if self.E // 2 > 2:
                return v
            v[self.lepton_flavor - 1] = self.lepton_number
        return v

    @property
    def is_neutrino(self):
        """return whether particle is a neutrino"""
        if self.is_lepton:
            f = self.E // 2
            if f == self.E_LEPTON_ANY - 1:
                return self.E - 2 * f == 0
            return self.Z == 0
        return False

    @property
    def lepton_flavor_vec_ext(self):
        """return lepton flavor vector [l_e, l_nu, l_tau, l_x]"""
        v = [0, 0, 0, 0]
        if self.is_lepton:
            if self.E // 2 == self.E_LEPTON_ANY - 1:
                return v
            v[self.lepton_flavor - 1] = self.lepton_number
        return v

    _isobar = re.compile(r"[Aa][=:-]?(\d+)")
    _isotone = re.compile(r"[Nn][=:-](\d+)")
    _element = re.compile(r"(?:[Zz][=:-]?|[eE])(\d+)")

    @classmethod
    def ion2bfzae_(cls, s="", element=False, isomer=False):
        pass

    @classmethod
    def ion2bfzae(cls, s="", element=False, isomer=False):
        """Interpret string and return component values."""
        if isinstance(s, (bytes, np.bytes_)):
            s = s.decode("us_ascii")
        s, sz, sa, se = cls._decompose(s, element=element)
        if sz == sa == se == "":
            return cls.VOID

        valid = True
        f = cls.F_ISOMER
        b = 0
        z = None
        try:
            a = int(sa)
        except ValueError:
            a = 0
            if sa != "":
                valid = False
        if sz in ("a=", "A=", "A:", "a:", "a-", "A-", "a", "A"):
            if a == 0:
                valid = False
            f = cls.F_ISOBAR
            se = 0
            z = 0
        if sz in ("n=", "N=", "N:", "n:", "n-", "N-"):
            if a == 0:
                valid = False
            f = cls.F_ISOTONE
            se = 0
            z = 0
        if sz in ("z=", "Z=", "Z:", "Z:", "z-", "Z-", "z", "Z", "e", "E"):
            valid = a > 0
            z = a
            a = 0
            try:
                z = int(sa)
            except ValueError:
                valid = False
            else:
                f = cls.F_ELEMENT
                se = 0
        if z is None and not (
            s.endswith((",", ";", "-", "=", ":"))
            or s.startswith((",", ";", "-", "=", ":"))
        ):
            try:
                z = Elements.index(sz.strip(",;-=:"))
            except ValueError:
                pass
        if z is None:
            try:
                z = elements.index(sz.lower())
            except ValueError:
                pass
        if z is None and not (s.startswith("-") or (s.endswith("-"))):
            try:
                z = z2Name.index(sz.strip("-"))
            except ValueError:
                pass
        if z is None and not (s.startswith("-") or (s.endswith("-"))):
            try:
                z = z2name.index(sz.strip("-"))
            except ValueError:
                z = 0
                if sz != "":
                    valid = False
        try:
            e = int(se)
        except ValueError:
            e = 0
            if not isomer:
                f = cls.F_ISOTOPE
            if se != "":
                valid = False
        if sz in ("nt", "Nt") and a == 0:
            f = cls.F_ELEMENT
        elif sz in ("n", "N") and a > 1:
            z = 7
        elif sz in ("p", "P") and a > 1:
            z = 15
        elif sz.lower() in ("g", "gam", "gamma") and a > 0:
            f = cls.F_BOSON
            e = a
            a = 0
            valid = True
        elif z == a == 0 and e > 0:
            f = cls.F_BOSON
        else:
            try:
                i = cls.SPECIAL[s]
                b, f, z, a, e = cls.idx2bfzae(i)
                valid = True
            except:
                pass
        if z is None:
            z = 0
            valid = False
        if f == cls.F_ISOTOPE and z > 0 and a == 0:
            f = cls.F_ELEMENT
        if f == cls.F_ISOMER and z > 0 and a == 0:
            valid = False
        if not valid:
            b, f, z, a, e = cls.VOID
        return b, f, z, a, e

    def isomer(self, Z=None, A=None, N=None, E=None):
        """
        return self if isomer, isomer if enough info is provided, or VOID otherwise

        """
        if self.is_isomer and (E is None or isinstance(E, Mapping)):
            if A is N is Z is None:
                return self
            E = self.E
        x = self.isotope(A=A, Z=Z, N=N)
        if x is VOID:
            return ion
        if E is None:
            E = isomermap
        if isinstance(E, Mapping):
            try:
                E = E[x]
            except KeyError:
                # should default behavior be to set E = 0 instead?
                return VOID
        if E is None:
            return VOID
        return Isomer(A=x.A, Z=x.Z, E=E)

    def isotope(self, Z=None, A=None, N=None):
        """
        return isotope if isomer, or self if isotope, or VOID otherwise
        """
        if self.is_isotope and Z == A == N is None:
            return self

        n = 3 - [Z, A, N].count(None)
        if self.is_isotope or self.is_isomer:
            if n < 2 and Z is None:
                Z = self.Z
                n += 1
            if n < 2 and A is None:
                A = self.A
                n += 1
        elif self.is_isobar:  # has A
            if n < 2 and A is None:
                A = self.A
        elif self.is_isotone:  # has N
            if n < 2 and N is None:
                N = self.N
        elif self.is_element:  # has Z
            if n < 2 and Z is None:
                Z = self.Z

        # create output
        if Z is None:
            if A is None:
                return VOID
            else:
                if N is None:
                    return VOID
                else:
                    if N > A:
                        return VOID
                    return Isotope(N=N, A=A)
        else:
            if A is None:
                if N is None:
                    return VOID
                else:
                    return Isotope(Z=Z, N=N)
            else:
                if N is None:
                    if Z > A:
                        return VOID
                    return Isotope(Z=Z, A=A)
                else:
                    if N + Z == A:
                        return Isotope(Z=Z, A=A)
        return VOID

    def isobar(self, A=None):
        """
        return isobar if isomer or isotope, or self if isobar, or VOID otherwise
        """
        if A is not None:
            return Isobar(A=A)
        if self.is_isobar:
            return self
        if self.is_isotope or self.is_isomer:
            return Isobar(A=self.A)
        return self.VOID

    def isotone(self, N=None):
        """
        return isotone if isomer or isotope, or self if isotone, or VOID otherwise
        """
        if N is not None:
            return Isotone(N=N)
        if self.is_isotone:
            return self
        if self.is_isotope or self.is_isomer:
            return Isotone(N=self.N)
        return self.VOID

    def element(self, Z=None):
        """
        return element if isomer or isotope, or self if element, or VOID otherwise
        """
        if Z is not None:
            return Element(Z=Z)
        if self.is_element:
            return self
        if self.is_isotope or self.is_isomer:
            return Element(Z=self.Z)
        return self.VOID

    _translate_nugrid = {
        "cd*15": "cd115m",
        "lu*76": "lu115m",
        "tag80": "ta180g",
        "prot": "pn1",
        "neut": "nt1",
        "OOOOO": "g",
    }

    _translate_limongi = {
        "alg6": "al26g",
        "alm6": "al26m",
    }

    _translate_reaclib = {
        "al-6": "al26g",
        "al*6": "al26m",
    }

    _translate_extra = {
        "alpha": "he4",
        "gamma": "g",
        "nueb": "nbe",
        "numb": "nbm",
        "nutb": "nbt",
        "e": "e-",
        "bp": "e+",
        "bm": "e-",
        "tau": "t-",
        "tau-": "t-",
        "tau+": "t+",
        "m": "m-",
        "mu": "m-",
        "mu-": "m-",
        "mu+": "m+",
    }

    _translate = {}
    _translate.update(_translate_nugrid)
    _translate.update(_translate_reaclib)
    _translate.update(_translate_limongi)
    _translate.update(_translate_extra)

    _html = re.compile(r"<sup>([0-9mg*]+)</sup>([a-zA-Z]+)")
    _kobayashi = re.compile(r"\^(\d+)\^([A-Za-z])")

    @classmethod
    def _decompose(
        cls,
        s="",
        element=False,
    ):
        """
        This routine TRIES to decompose string in
        1) element symbol
        2) mass number
        3) excited state number
        All are returned as strings.

        The optional parameter 'element' enforces that
        'p' is taken as *element* phosphorus (not H1)
        'n' is taken as *element* nitrogen (not neutron)
        """

        s = s.strip()

        x = cls._html.findall(s)
        if len(x) > 0:
            s = "".join(x[0][::-1])
        x = cls._kobayashi.findall(s)
        if len(x) > 0:
            s = "".join(x[0][::-1])

        s = cls._translate.get(s.lower(), s)

        name = s.strip()
        el = ""
        a = ""
        e = ""

        # get numbers
        n = re.findall(r"\d+", name)

        # get strings
        cxx = re.findall(r"\D+", name)

        c = list()
        for x in cxx:
            xx = x.split("-")
            cy = [y for y in xx if y != ""]
            c += cy
        cx = c
        c = list()
        for cz in cx:
            cy = ""
            while cz.startswith(iso_string_eany):
                cy += cz[0]
                cz = cz[1:]
            if len(cy) > 0:
                c.append(cy)
            if len(cz) > 0:
                c.append(cz)
        cx = c
        c = list()
        for cz in cx:
            cy = ""
            while cz.endswith(iso_string_eany):
                cy += cz[-1]
                cz = cz[:-1]
            if len(cz) > 0:
                c.append(cz)
            if len(cy) > 0:
                c.append(cy)

        if len(c) == 1 and len(c[0]) > 1:
            cx = str.lower(c[0])
            if cx[0] in iso_strings_:
                if cx not in elements and cx[1:] in elements_:
                    c = [c[0][0], c[0][1:]]

            if len(c) == 1 and cx[-1] in iso_strings_:
                if cx not in elements and cx[:-1] in elements_:
                    c = [c[0][:-1], c[0][-1]]

            if len(c) == 1 and cx[-1] in iso_strings_:
                if cx in elements and cx[:-1] in elements_:
                    if 0 < len(n) <= 2:
                        nx = int(n[0])
                        zx = el2z[cx]
                        if not cls._valid_isotope(nx, zx):
                            c = [c[0][:-1], c[0][-1]]

            if (
                len(c) == 1
                and cx[0] in iso_strings_
                and not cx[1] in iso_special_particle
            ):
                if cx in elements and cx[1:] in elements_:
                    if len(n) == 1:
                        nx = int(n[0])
                        zx = el2z[cx]
                        if not cls._valid_isotope(nx, zx):
                            c = [c[0][1:], c[0][0]]

        if len(c) == 2:
            if c[0] in iso_strings_:
                c = c[::-1]
            if c[0][0] == iso_string_eany:
                if len(n) > 1:
                    c[0] = "m" + c[0][1:]
                c = c[::-1]
            if c[0] in elements_ and not (
                c[1] in iso_strings or c[1] == iso_string_eany
            ):
                return ("-",) + ("",) * 3

        # this filter may come later
        if len(c) > 2 and len(n) > 1:
            return ("-",) + ("",) * 3

        if len(n) > 0:
            a = n[0]
        if len(n) > 1:
            e = n[1]
        if len(n) > 2:
            return ("-",) + ("",) * 3
            # raise ValueError(f"Can't understand isotope '{s}'.")
        if len(c) > 0:
            el = c[0]
        if len(c) == 2 and c == ["(", ")"]:
            if len(n) == 1:
                a = n[0]
                el = "Z="
                e = ""
                c = []
                n = []
            else:
                return (s,) + ("",) * 3

        # experimental
        if el in iso_special_mass and e == "":
            e = a
            a = str(iso_special_mass[el])
            if el in iso_special_element and len(n) > 0:
                if n[0] != a:
                    return ("-",) + ("",) * 3
            elif el in iso_special_particle:
                if len(n) == 1:
                    if len(c) == 1:
                        return ("-",) + ("",) * 3
                elif len(n) > 1:
                    return ("-",) + ("",) * 3
            n.insert(0, a)
            c[0] = el = elements[el2z[el]]

        if len(c) == 2:
            if c[1] in iso_strings_ground:
                e = "0"
                if len(n) > 1:
                    return (s,) + ("",) * 3
            elif c[1] in iso_strings_excite and len(n) == 1:
                e = "1"
            elif c[1][0] == iso_string_eany and len(n) == 1:
                e = str(len(c[1]))
                assert c[1].count(iso_string_eany) == len(c[1])
                if e == "1":
                    e = str(cls.EANY)
            if not c[1] in iso_strings and not c[1][0] == iso_string_eany:
                return (s,) + ("",) * 3

        if len(c) == 1 and c[0][-1] == iso_string_eany:
            e = 0
            while c[0][-1] == iso_string_eany:
                c[0] = c[0][:-1]
                e += 1
            assert e == 1
            e = str(e)
            el = c[0]

        if len(c) == 1 and c[0][0] == iso_string_eany:
            e = 0
            while c[0][0] == iso_string_eany:
                c[0] = c[0][1:]
                e += 1
            assert e == 1
            e = str(e)
            el = c[0]

        if s == "a" and a == "":
            el = "He"
            a = "4"
        # this is a possible conflict with potassium
        elif (element) and s == "p":
            el = "P"
        elif s == "p":
            el = "H"
            a = "1"
        elif el in ("p", "pn") and a == "1":
            el = "H"
        elif s == "pn":
            el = "H"
        elif el in ("d", "D"):
            el = "H"
            if a not in ("", "2"):
                return ("-",) + ("",) * 3
                # raise AttributeError('"d" already implies mass; if supplied needs to be "2".')
            a = "2"
        elif el in ("t", "T"):
            el = "H"
            if a not in ("", "3"):
                return ("-",) + ("",) * 3
                # raise AttributeError('"t" already implies mass; if supplied needs to be "3"')
            a = "3"
        elif (element) and s == "n":
            el = "N"
        elif s == "n":
            el = "nt"
            a = "1"
        elif el in ("n", "nt") and a == "1":
            el = "nt"
        elif s in ("g", "G"):
            el = ""
            a = ""
            e = "1"
        elif s == "b":
            s = el = "e-"
        elif s.lower() in ("e-", "b-", "bd", "pc"):
            s = el = "e-"
        elif (s.lower() in ("e+", "b+", "ec")) or (
            (not element) and (s in ("pd", "PD"))
        ):
            s = el = "e+"
        elif (not element) and (s.lower() == "ps"):
            s = "h1"
            a = "1"
            el = "h"
        elif (not element) and (s.lower() == "ns"):
            s = "nt1"
            a = "1"
            el = "nt"
        el = el.strip()
        #        if len(el) == 2 and el(2)
        a = a.strip()
        e = e.strip()

        if s.endswith("-") and not el.endswith("-"):
            el = s

        return s, el, a, e

    def _deexcite(self):
        """Do not use."""
        self.E = 0
        self.idx = self.bfzae2idx(F=self.F, Z=self.Z, A=self.A, E=self.E, B=self.B)

    @property
    def deexcited(self):
        return ion(
            F=self.F,
            A=self.A,
            Z=self.Z,
            B=self.B,
        )

    def __str__(self):
        return self._name(upcase=True)

    def __format__(self, format_spec):
        return "{:{}}".format(self.__str__(), format_spec)

    def nugridpy(self):
        # need to update for isomeres
        return f"{Elements[self.Z]}-{self.A:d}"

    _translate_to_mesa = {
        "n": "neut",
    }

    @property
    def mesa(self):
        s = self.name()
        return self._translate_to_mesa.get(s, s)

    @property
    def kobayashi(self):
        # need to update for isomeres
        return f"^{self.A:d}^{Elements[self.Z]}"

    _custom_add = False

    def _add(self, x, sign1=+1, sign2=+1):
        if np.shape(x) != ():
            return NotImplemented
        if not self._is_ion(x):
            y = ion(x)
            if y.idx == self.VOID_IDX:
                raise AttributeError(
                    f"Cannot convert {x!r} to {Ion.__class__.__name__!r}."
                )
            else:
                x = y
        if x._custom_add:
            return NotImplemented
        A = sign1 * self.A + sign2 * x.A
        Z = sign1 * self.Z + sign2 * x.Z
        E = sign1 * self.E + sign2 * x.E
        N = sign1 * self.N + sign2 * x.N

        if not (
            ((self.F & self.F_GROUP_MASK == self.F_ISOMER) and x.is_photon)
            or (x.F & self.F_GROUP_MASK == self.F_ISOMER)
            and x.is_photon
        ):
            E = None
        if (self.is_lepton and not x.is_lepton) or (x.is_lepton and not self.is_lepton):
            if self.F_GROUP_MASK == self.F_ISOMER or x.F_GROUP_MASK == self.F_ISOMER:
                E = 0
            else:
                E = None
            N = A - Z
        elif (self.F & self.F_GROUP_MASK == self.F_ELEMENT) and (
            x.F & self.F_GROUP_MASK == self.F_ELEMENT
        ):
            A = None
            N = None
        elif (self.F & self.F_GROUP_MASK == self.F_ISOBAR) and (
            x.F & self.F_GROUP_MASK == self.F_ISOBAR
        ):
            Z = None
            N = None
        elif (self.F & self.F_GROUP_MASK == self.F_ISOTONE) and (
            x.F & self.F_GROUP_MASK == self.F_ISOTONE
        ):
            Z = None
            A = None
        elif (self.F & self.F_GROUP_MASK == self.F_ISOBAR) or (
            x.F & self.F_GROUP_MASK == self.F_ISOBAR
        ):
            N = A - Z
        else:
            A = N + Z
        return ion(
            Z=Z,
            A=A,
            E=E,
            N=N,
        )

    def __add__(self, x):
        return self._add(x)

    __radd__ = __add__

    def __sub__(self, x):
        return self._add(x, +1, -1)

    def __rsub__(self, x):
        return self._add(x, -1, +1)

    def __mul__(self, x):
        if np.shape(x) != ():
            return NotImplemented
        return self.__class__(
            Z=self.Z * x,
            A=self.A * x,
            E=self.E,
            F=self.F,
        )

    __rmul__ = __mul__

    def __neg__(self):
        return self.__class__(
            Z=-self.Z,
            A=-self.A,
            E=self.E,
            F=self.F,
            B=self.B,
        )

    def __bool__(self):
        return self.idx != self.VOID_IDX

    def __hash__(self):
        return self.idx

    def __lt__(self, x):
        if np.shape(x) != ():
            return NotImplemented
        if not isinstance(x, Ion):
            x = self.factory(x)
        return self.idx < x.idx

    def __eq__(self, x):
        if np.shape(x) != ():
            return NotImplemented
        if not isinstance(x, Ion):
            x = self.factory(x)
        return self.idx == x.idx

    def __le__(self, x):
        if np.shape(x) != ():
            return NotImplemented
        if not isinstance(x, Ion):
            x = self.factory(x)
        return self.idx <= x.idx

    def __gt__(self, x):
        return not self.__le__(x)

    def __ge__(self, x):
        return not self.__lt__(x)

    def __len__(self):
        return 1

    def __repr__(self):
        base = f"('{self.Name():s}')"
        return self.__class__.__name__ + base

    def _repr_latex_(self):
        s = self.LaTeX(math=False)
        return self.__class__.__name__ + "(" + s + ")"

    def __getstate__(self):
        # need to add version number
        # should save index instead
        #        print('xxx_save')
        return (
            self.B,
            self.F,
            self.Z,
            self.A,
            self.E,
            self.N,
            self.idx,
        )

    def __setstate__(self, x):
        # need to add check of version number
        #        print('xxx_load')
        self.__init__(_load=x)

    @property
    def element_name(self):
        return z2name[self.Z]

    @property
    def is_nucleus(self):
        return self.is_isomer or self.is_isotope

    @property
    def is_lepton(self):
        return self.F & self.F_GROUP_MASK == self.F_LEPTON and self.F != -1

    @property
    def is_isomer(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOMER

    @property
    def is_isotope(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTOPE

    @property
    def is_isobar(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOBAR

    @property
    def is_isotone(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTONE

    @property
    def is_element(self):
        return self.F & self.F_GROUP_MASK == self.F_ELEMENT

    @property
    def is_ionize(self):
        return self.F & self.F_GROUP_MASK == self.F_IONIZE

    @property
    def is_boson(self):
        return self.F & self.F_GROUP_MASK == self.F_BOSON

    @property
    def is_photon(self):
        return (
            (self.F & self.F_GROUP_MASK == self.F_BOSON)
            and (self.Z == 0)
            and (self.A == 0)
            and (self.E > 0)
        )

    @property
    def is_void(self):
        return self.F == -1

    @property
    def is_valid(self):
        return self.F != -1

    @property
    def is_vac(self):
        return self.F == 0

    @property
    def type(self):
        if self.is_lepton:
            return "lepton"
        if self.is_isotope:
            return "isotope"
        if self.is_isobar:
            return "isobar"
        if self.is_isotone:
            return "isotone"
        if self.is_element:
            return "element"
        if self.is_isomer:
            return "isomer"
        if self.is_ionize:
            return "ionize"
        if self.is_nucleon:
            return "nucleon"
        if self.is_photon:
            return "photon"
        if self.is_boson:
            return "boson"
        if self.is_void:
            return "void"
        if self.is_vac:
            return "vac"
        return "unknown"

    @property
    def anti_particle(self):
        return self.__class__(
            A=self.A,
            Z=-self.Z,
            F=self.F,
            B=self.B,
        )

    def R(self, I=None, O=None):
        """
        Reaction rate function.
        If only one parameter is given it is interpreted as 'outgoing'
        If no parameter is provided, just a deep copy of self is returned
        """
        if O is None:
            if I is None:
                return self
            return self(None, I)
        return self(I, O)

    def __call__(self, *args):
        """
        implement reaction semantics

        if only one parameter is supplied, assume it is the output
        channel

        if more than one is supplied, assume the first is the in
        channel, the other are the out channel

        in/out channels may be iterables

        Examples
        ion('c12')('p','g') --> Isotope('N13')

        ion('c12')('p') -->


        TODO - add tuples to group channels
        """
        result = self
        if len(args) == 1:
            if np.isscalar(args[0]) or args[0] is None:
                result -= args[0]
            else:
                for i in args[0]:
                    result -= i
            return result
        if np.isscalar(args[0]) or args[0] is None:
            result += args[0]
        else:
            for i in args[0]:
                result += i
        for i in args[1:]:
            if np.isscalar(i) or i is None:
                result -= i
            else:
                for j in i:
                    result -= j
        return result


ufunc_A = lambda y: np.array(np.frompyfunc(lambda x: x.A, 1, 1)(y), dtype=int)
ufunc_Z = lambda y: np.array(np.frompyfunc(lambda x: x.Z, 1, 1)(y), dtype=int)
ufunc_N = lambda y: np.array(np.frompyfunc(lambda x: x.N, 1, 1)(y), dtype=int)
ufunc_E = lambda y: np.array(np.frompyfunc(lambda x: x.E, 1, 1)(y), dtype=int)
ufunc_F = lambda y: np.array(np.frompyfunc(lambda x: x.F, 1, 1)(y), dtype=int)
ufunc_B = np.frompyfunc(lambda x: x.B, 1, 1)  # these may be larger than int64

ufunc_isomer = lambda y: np.array(
    np.frompyfunc(lambda x: x.isomer(), 1, 1)(y), dtype=object
)
ufunc_isotope = lambda y: np.array(
    np.frompyfunc(lambda x: x.isotope(), 1, 1)(y), dtype=object
)
ufunc_element = lambda y: np.array(
    np.frompyfunc(lambda x: x.element(), 1, 1)(y), dtype=object
)
ufunc_isobar = lambda y: np.array(
    np.frompyfunc(lambda x: x.isobar(), 1, 1)(y), dtype=object
)
ufunc_isotone = lambda y: np.array(
    np.frompyfunc(lambda x: x.isotone(), 1, 1)(y), dtype=object
)

ufunc_element_name = lambda y: np.array(
    np.frompyfunc(lambda x: x.element_name, 1, 1)(y)
)
ufunc_element_symbol = lambda y: np.array(
    np.frompyfunc(lambda x: x.element_symbol(upase=True), 1, 1)(y)
)
ufunc_element_sym_lc = lambda y: np.array(
    np.frompyfunc(lambda x: x.element_symbol(upase=False), 1, 1)(y)
)

ufunc_isomer_idx = lambda y: np.array(
    np.frompyfunc(lambda x: x.isomer().idx, 1, 1)(y), dtype=int
)
ufunc_isotope_idx = lambda y: np.array(
    np.frompyfunc(lambda x: x.isotope().idx, 1, 1)(y), dtype=int
)
ufunc_element_idx = lambda y: np.array(
    np.frompyfunc(lambda x: x.element().idx, 1, 1)(y), dtype=int
)
ufunc_isobar_idx = lambda y: np.array(
    np.frompyfunc(lambda x: x.isobar().idx, 1, 1)(y), dtype=int
)
ufunc_isotone_idx = lambda y: np.array(
    np.frompyfunc(lambda x: x.isotone().idx, 1, 1)(y), dtype=int
)

ufunc_is_ion = lambda y: np.array(
    np.frompyfunc(lambda x: isinstance(x, Ion), 1, 1)(y), dtype=bool
)

ufunc_is_lepton = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_lepton, 1, 1)(y), dtype=bool
)
ufunc_is_isotope = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_isotope, 1, 1)(y), dtype=bool
)
ufunc_is_isobar = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_isobar, 1, 1)(y), dtype=bool
)
ufunc_is_isotone = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_isotone, 1, 1)(y), dtype=bool
)
ufunc_is_element = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_element, 1, 1)(y), dtype=bool
)
ufunc_is_isomer = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_isomer, 1, 1)(y), dtype=bool
)
ufunc_is_ionize = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_ionize, 1, 1)(y), dtype=bool
)
ufunc_is_nucleus = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_nucleus, 1, 1)(y), dtype=bool
)
ufunc_is_photon = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_photon, 1, 1)(y), dtype=bool
)
ufunc_is_void = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_void, 1, 1)(y), dtype=bool
)
ufunc_is_vac = lambda y: np.array(
    np.frompyfunc(lambda x: x.is_vac, 1, 1)(y), dtype=bool
)

ufunc_type = lambda y: np.array(np.frompyfunc(lambda x: x.type, 1, 1)(y), dtype=np.str)

ufunc_idx = lambda y: np.array(np.frompyfunc(lambda x: x.idx, 1, 1)(y), dtype=int)
ufunc_idx_ZA = lambda y: ((y // Ion.ZMUL) % Ion.ZMAX, (y // Ion.AMUL) % Ion.AMAX)
_ufunc_idx_ZAE = lambda y: (
    (y // Ion.ZMUL) % Ion.ZMAX,
    (y // Ion.AMUL) % Ion.AMAX,
    (y // Ion.EMUL) % Ion.EMAX,
)
_ufunc_idx_ZAEF = lambda y: (
    (y // Ion.ZMUL) % Ion.ZMAX,
    (y // Ion.AMUL) % Ion.AMAX,
    (y // Ion.EMUL) % Ion.EMAX,
    (y // Ion.FMUL) % Ion.FMAX,
)


def ufunc_idx_ZN(idx):
    Z, A = ufunc_idx_ZA(idx)
    return Z, A - Z


def ufunc_idx_ZAE(idx):
    Z, A, E, F = _ufunc_idx_ZAEF(idx)
    ii = F == Ion.F_ISOMER
    E[~ii] = -1
    return Z, A, E


def ufunc_idx_ZNE(idx):
    Z, A, E, F = _ufunc_idx_ZAEF(idx)
    ii = F == Ion.F_ISOMER
    E[~ii] = -1
    return Z, A - Z, E


ufunc_idx_from_ZA = (
    lambda Z, A: np.array(Z) * Ion.ZMUL
    + np.array(A) * Ion.AMUL
    + Ion.F_ISOTOPE * Ion.FMUL
)


ufunc_idx_from_ZN = (
    lambda Z, N: np.array(Z) * Ion.ZMUL
    + (np.array(N) + np.array(N)) * Ion.AMUL
    + Ion.F_ISOTOPE * Ion.FMUL
)


def ufunc_idx_isomer_from_ZAE(Z, A, E):
    E = np.array(E)
    F = np.full(Z.shape, Ion.F_ISOMER)
    ii = E == -1
    F[ii] = Ion.F_ISOTOPE
    E[ii] = 0
    return np.array(Z) * Ion.ZMUL + np.array(A) * Ion.AMUL + F * Ion.FMUL


def ufunc_idx_isomer_from_ZNE(Z, N, E):
    Z = np.array(Z)
    N = np.array(N)
    E = np.array(E)
    F = np.full(Z.shape, Ion.F_ISOMER)
    ii = E == -1
    F[ii] = Ion.F_ISOTOPE
    E[ii] = 0
    return Z * Ion.ZMUL + (N + Z) * Ion.AMUL + F * Ion.FMUL


ufunc_ion = lambda y: np.array(np.frompyfunc(lambda x: ion(x), 1, 1)(y), dtype=object)


ufunc_ion_from_ZA = lambda Z, A: np.array(
    np.frompyfunc(lambda z, a: ion(Z=z, A=a), 2, 1)(Z, A), dtype=object
)


ufunc_isotope_from_ZA = lambda Z, A: np.array(
    np.frompyfunc(lambda z, a: ion(Z=z, A=a, isotope=True), 2, 1)(Z, A), dtype=object
)


ufunc_ion_from_ZN = lambda Z, N: np.array(
    np.frompyfunc(lambda z, n: ion(Z=z, N=n), 2, 1)(Z, N), dtype=object
)


ufunc_isotope_from_ZN = lambda Z, N: np.array(
    np.frompyfunc(lambda z, n: ion(Z=z, N=n, isotope=True), 2, 1)(Z, N), dtype=object
)


ufunc_ion_from_ZAE = lambda Z, A, E: np.array(
    np.frompyfunc(lambda z, a, e: ion(Z=z, A=a, E=e), 3, 1)(Z, A, E),
    dtype=object,
)


ufunc_isomer_from_ZAE = lambda Z, A, E: np.array(
    np.frompyfunc(lambda z, a, e: ion(Z=z, A=a, E=e, isomer=True), 3, 1)(Z, A, E),
    dtype=object,
)


ufunc_ion_from_ZNE = lambda Z, N, E: np.array(
    np.frompyfunc(lambda z, n, e: ion(Z=z, N=n, E=e), 3, 1)(Z, N, E),
    dtype=object,
)


ufunc_isomer_from_ZNE = lambda Z, N, E: np.array(
    np.frompyfunc(lambda z, n, e: ion(Z=z, N=n, E=e, isomer=True), 3, 1)(Z, N, E),
    dtype=object,
)


ufunc_ion_from_kw = lambda y, kw: np.array(
    np.frompyfunc(lambda x: ion(x, **kw), 1, 1)(y), dtype=object
)


ufunc_ion_from_idx = lambda y: np.array(
    np.frompyfunc(lambda x: ion(idx=x), 1, 1)(y), dtype=object
)


# some convenience functions
def ionarr(arr):
    if isinstance(arr, str):
        arr = re.split("[^-0-9a-zA-Z*]+", arr)
    if not isinstance(arr, Iterable):
        arr = (arr,)
    return np.array([ion(a) for a in arr], dtype=object)


def ision(x):
    return isinstance(x, Ion)


def isionclass(x):
    return issubclass(x, Ion)


def get_ufunc(mode):
    return sys.modules[__name__].__dict__[f"ufunc_is_{mode}"]


def get_ufunc_idx(mode):
    return sys.modules[__name__].__dict__[f"ufunc_{mode}_idx"]


# the following classes need to be fleshed out and used consequnetly
# in particular, maybe some filter of constructor can be done?
class Element(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs["element"] = True
        super().__init__(*args, **kwargs)
        assert self.is_element


class Isotone(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs["isotone"] = True
        super().__init__(*args, **kwargs)
        assert self.is_isotone


class Isobar(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs["isobar"] = True
        super().__init__(*args, **kwargs)
        assert self.is_isobar


class Isotope(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_isotope


class Isomer(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs["isomer"] = True
        super().__init__(*args, **kwargs)
        assert self.is_isomer


class Photon(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_photon


class Lepton(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_lepton

    flavor = Ion.lepton_flavor
    flavor_vec = Ion.lepton_flavor_vec
    flavor_vec_ext = Ion.lepton_flavor_vec_ext
    number = Ion.lepton_number


class Ionize(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_ionize


class Void(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_void


class Vac(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_vac


# creation functions/factories
def element(Z):
    return Element(Z=Z)


def isotone(N):
    return Isotone(N=N)


def isobar(A):
    return Isobar(A=A)


def isotope(A, Z):
    return Isotope(A=A, Z=Z)


def isomer(A, Z, E=0):
    return Isomer(A=A, Z=Z, E=E)


def photon(E=1):
    return Photon(E=E)


def lepton(Z=0, number=1, flavor=None):
    """
    use with E=1 for lepton, E=0 for anti-lepton
    """
    if number is None:
        if Z == 1:
            number = -1
        elif Z == -1:
            number = 1
    if flavor is None:
        flavor = 1
    if flavor == 0:
        assert Z == 0
        E = Ion.E_LEPTON_FLAVOR_MUL * (Ion.E_LEPTON_ANY - 1) + number
    else:
        E = (number + 1) // 2
        assert 0 < flavor <= Ion.E_LEPTON_X
        assert Z == 0 or Z == 1 or (Z == -1 and E == 1)
        if Z == 1 and E == 1:
            Z = -1
        E += Ion.E_LEPTON_FLAVOR_MUL * (flavor - 1)
    return Lepton(Z=Z, E=E, lepton=True)


def ionize(Z=0, A=0, E=0):
    assert E <= Z, "too much ionization"
    assert (A == 0) or (A >= Z), "no such isotope"
    return Ionize(Z=Z, E=E, A=A, F=Ion.F_IONIZE)


def void():
    return Void(Ion.VOID)


class IonFactory:
    """
    Factory function class.

    This one should be used to allow for later refactoring in derived
    classes for different ion types.
    """

    _parms = (
        "Z",
        "A",
        "N",
        "E",
    )

    def __call__(self, *args, **kwargs):
        if len(args) == 0 and len(kwargs) == 0:
            return
        if len(args) == 1 and len(kwargs) == 0 and ision(args[0]):
            return args[0]
        if len(args) == 1 and isinstance(args[0], (np.ndarray, list, set)):
            return ufunc_ion_from_kw(args[0], kwargs)
        if len(args) == 1 and isinstance(args[0], dict):
            kwargs.update(args[0])
            args = ()

        # allow for ion(A=[...], Z=[...]) format
        parms = np.array([kwargs.get(k, None) for k in self._parms], dtype=object)
        if np.any(parms is not None):
            b = np.broadcast(*parms)
            if b.size > 1:
                ions = np.empty(b.shape, dtype=object)
                pl = list()
                vl = list()
                for k, v in zip(self._parms, parms):
                    if v is not None:
                        vl.append(np.broadcast_to(v, b.shape))
                        pl.append(k)
                for i, vx in enumerate(zip(*(vf.flat for vf in vl))):
                    kw = kwargs.copy()
                    for p, v in zip(pl, vx):
                        kw[p] = v
                    ions.flat[i] = ion(*args, **kw)
                return ions

        kwargs["factory"] = True
        ix = Ion(*args, **kwargs)
        if ix.is_isomer:
            return Isomer(ix)
        if ix.is_isotope:
            return Isotope(ix)
        if ix.is_element:
            return Element(ix)
        if ix.is_isobar:
            return Isobar(ix)
        if ix.is_isotone:
            return Isotone(ix)
        if ix.is_photon:
            return Photon(ix)
        if ix.is_lepton:
            return Lepton(ix)
        if ix.is_ionize:
            return Ionize(ix)
        if ix.is_void:
            return Void(ix)
        if ix.is_vac:
            return Vac(ix)
        debug = kwargs.get("debug", True)
        if debug:
            raise Exception("Unknown Ion type.")
        return ix

    def __getattr__(self, attr):
        if not isinstance(attr, Ion):
            x = ion(attr)
        else:
            x = attr
        if x != Ion.VOID:
            return self(x)
        raise AttributeError(attr)


ion = IonFactory()
Ion.factory = ion

VAC = ion(Ion.VAC)
GAMMA = ion(Ion.GAMMA)
VOID = ion(Ion.VOID)
NEUTRON = ion(Z=0, A=1)
PROTON = ion(Z=1, A=1)
ELECTRON = lepton(Z=-1, number=1, flavor=1)
POSITRON = lepton(Z=+1, number=-1, flavor=1)
NUE = lepton(Z=0, number=1, flavor=1)
NBE = lepton(Z=0, number=-1, flavor=1)
NU = lepton(Z=0, number=0, flavor=0)
ALPHA = ion(Z=2, A=4)
VOID_IDX = Ion.VOID_IDX


class IsomerMap(collections.defaultdict):
    """
    define map that can be used as 'E' parameter in ion.isomer function

    TODO - add more isotopes, provide function to load from file
    """

    def __init__(self, default=lambda: 0, isomap=None):
        super().__init__(default)
        if isomap is not None:
            self.update(isomap)

    def __getitem__(self, item):
        if not isinstance(item, Ion):
            item = ion(item)
        return super().__getitem__(item)

    def __setitem__(self, item, value):
        if not isinstance(item, Ion):
            item = ion(item)
        super().__setitem__(item, value)


isomermap = IsomerMap(isomap={Isotope(A=180, Z=73): 1})

# registry of other bits [set, but not yet used except to check for duplicates]
other_bits_register = {}


def register_other_bits(bit, name=None, class_=None, overwrite=True):
    if isinstance(bit, type):
        class_ = bit
        name = class_.__name__
        bit = class_.__dict__["F_OTHER_" + name.upper()]

    if name is not None and class_ is None:
        name, class_ = class_, name
    if name is None:
        name = class_.__name__

    # this may need to become more sophisticated to allow re-compilation
    if not overwrite:
        assert bit not in other_bits_register, "bit already exists"

    other_bits_register[bit] = (name, class_)


def __getattr__(name):
    if name.startswith("__"):
        m = sys.modules[__name__].__dict__
        if name in m:
            return m[name]
        else:
            raise AttributeError()
    if name == "ioncacheza":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheZA(clone=m.get("ioncachezn", None))
        return m[name]
    elif name == "ioncachezn":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheZN(clone=m.get("ioncacheza", None))
        return m[name]
    elif name == "ioncachezae":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheZAE(clone=m.get("ioncachezne", None))
        return m[name]
    elif name == "ioncachezne":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheZNE(clone=m.get("ioncachezae", None))
        return m[name]
    elif name == "ioncachename":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheName(clone=m.get("ioncachename", None))
        return m[name]
    elif name == "ioncacheidx":
        m = sys.modules[__name__].__dict__
        if name not in m:
            m[name] = IonCacheIdx()
        return m[name]
    try:
        i = ion(name)
    except:
        raise AttributeError()
    if i is VOID:
        raise AttributeError()
    return i


# TODO - move into submodule, add load/share option
# TODO - add stringcache (dictionary)


class IonCacheIdx(object):
    def __init__(self, clone=None):
        if clone is None:
            self.db = dict()
        else:
            self.db = clone.copy()

    def __call__(self, idx):
        ions = np.ndarray(np.shape(idx), dtype=object)
        for i, ix in enumerate(np.array(idx).flat):
            io = self.db.get(ix, None)
            if io is None:
                io = ion(idx=ix)
                self.db[ix] = io
            ions.flat[i] = io
        return ions


class IonCacheZAE(object):
    def __init__(self, zmax=None, nmax=None, emax=None, clone=None):
        if clone is None:
            if zmax is None:
                zmax = 138
            if nmax is None:
                nmax = 256
            if emax is None:
                emax = 4
            self.ioncache0 = np.full(
                ((nmax + 1) * (zmax + 1) * (emax + 1)), False, dtype=bool
            )
            self.ioncache1 = np.ndarray(self.ioncache0.shape, dtype=object)
        else:
            if zmax is None:
                zmax = clone.zmax
            else:
                assert self.zmax == clone.zmax
            if nmax is None:
                nmax = clone.nmax
            else:
                assert self.nmax == clone.nmax
            if emax is None:
                emax = clone.emax
            else:
                assert self.emax == clone.emax
            self.ioncache0 = clone.ioncache0
            self.ioncache1 = clone.ioncache1
        self.zmax = zmax
        self.nmax = nmax
        self.emax = emax
        self._c1 = (self.nmax - 1) * self.emax
        self._c2 = self.emax
        self._c3 = self.nmax * self.emax + 1
        self._ufunc = ufunc_ion_from_ZAE

    def __call__(self, iz, ia, ie):
        s = np.shape(iz)
        iz, ia, ie = np.atleast_1d(iz, ia, ie)
        ix = self._c1 * iz + self._c2 * ia + ie + self._c3
        i0 = ~self.ioncache0[ix]
        (ii,) = np.where(i0)
        if len(ii) > 0:
            ions = self._ufunc(iz[ii], ia[ii], ie[ii])
            self.ioncache1[ix[ii]] = ions
            self.ioncache0[ix[ii]] = True
        return self.ioncache1[ix].reshape(s)[()]


class IonCacheZNE(IonCacheZAE):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._c1 = self.nmax * self.emax
        self._c2 = self.emax
        self._c3 = self._c1 + self._c2 + 1
        self._ufunc = ufunc_ion_from_ZNE


class IonCacheZA(object):
    def __init__(self, zmax=None, nmax=None, clone=None):
        if clone is None:
            if zmax is None:
                zmax = 138
            if nmax is None:
                nmax = 256
            self.ioncache0 = np.full(nmax * zmax, False, dtype=bool)
            self.ioncache1 = np.ndarray(self.ioncache0.shape, dtype=object)
        else:
            if zmax is None:
                zmax = clone.zmax
            else:
                assert self.zmax == clone.zmax
            if nmax is None:
                nmax = clone.nmax
            else:
                assert self.nmax == clone.nmax
            self.ioncache0 = clone.ioncache0
            self.ioncache1 = clone.ioncache1
        self.zmax = zmax
        self.nmax = nmax
        self._c1 = nmax - 1
        self._c3 = nmax
        self._ufunc = ufunc_ion_from_ZA

    def __call__(self, iz, ia, ie=-1):
        s = np.shape(iz)
        iz, ia, ie = np.atleast_1d(iz, ia, ie)
        if np.any(ie != -1):
            raise AttributeError("{ie=} not supported")
        ix = self._c1 * iz + ia + self._c3
        i0 = ~self.ioncache0[ix]
        (ii,) = np.where(i0)
        if len(ii) > 0:
            ions = self._ufunc(iz[ii], ia[ii])
            self.ioncache1[ix[ii]] = ions
            self.ioncache0[ix[ii]] = True
        return self.ioncache1[ix].reshape(s)[()]


class IonCacheZN(IonCacheZA):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._c1 = self.nmax
        self._c3 = self._c1 + 1
        self._ufunc = ufunc_ion_from_ZN


class IonCacheName(object):
    def __init__(self, clone=None):
        if clone is not None:
            self.ioncache = clone.ioncache.copy()
        else:
            self.cache = dict()
        self._ufunc = ufunc_ion

    def __call__(self, name):
        name = np.asarray(name)
        ions = np.ndarray(name.shape, dtype=object)
        for j, n in enumerate(name.flat):
            i = self.cache.get(n, None)
            if i is None:
                i = ion(n)
                self.cache[n] = i
            ions.flat[j] = i
        return ions[()]
