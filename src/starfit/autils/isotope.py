"""
Collection of routines for isotopes, most prominently the Ion class.
"""

# contains all the isotope routines

import collections
import re
from collections.abc import Iterable, Mapping

import numpy as np

from .utils import CachedAttribute, MetaSingletonHash

# Would b better if this was a class so we could write
# elements.H
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
    "ut",
    "fl",
    "up",
    "lv",
    "us",
    "uo",
)
# should be Uus Uuo but I will not use 3 letters
# as this breaks too many things

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
    "Flerovium",
    "(Ununquadium)",
    "livermorium",
    "(Ununhexium)",
    "(Ununseptium)",
    "(Ununoctium)",
)
z2name = tuple((name.lower() for name in z2Name))


# el2z provides easy convesion from element symbol to charge number
el2z = {el: i for i, el in enumerate(elements)}
el2z.update({el: i for i, el in enumerate(Elements[1:], start=1)})


def el2name(sym):
    return z2name[el2z[sym]]


class Ion(object, metaclass=MetaSingletonHash):
    """
    This is a generic class for 'ions'.
    It is designed to handel elements, isotioes, isomers, and isobars.
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

    g    = ( 0, 0, 0, 0)
    e-   = ( 1,-1, 0, 0)
    e+   = ( 1,+1, 0, 0)
    void = (-1, 0, 0, 0)

    We always put E1 as 'm'
    We destinguish isotopes and isomers - Al26g is not Al26
    """

    # in the future we want to use these
    # should flags be first or last?
    # maybe last?
    AMAX = 512
    ZMAX = 128
    EMAX = 256
    FMAX = 128  # flags - 7 bit
    # leave the sign bit for int64

    EMUL = 1
    AMUL = EMAX
    ZMUL = AMUL * AMAX
    FMUL = ZMUL * ZMAX
    BMUL = FMUL * FMAX
    # what lies Beyond ... (extra information for subclasses)
    # if these are used, the index no longer fits into int64

    # bit 0-2
    F_BOSON = 0
    F_ISOMER = 1
    F_ISOTOPE = 2
    F_ELEMENT = 3
    F_ISOBAR = 4
    F_ISOTONE = 5
    F_HADRON = 6  # currently not used
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
    VOID = (0, -1, 0, 0, 0)
    VOID_IDX = -1 * FMUL
    VOID_STRING = "-"
    GAMMA = (0, F_BOSON, 0, 0, +1)

    SPECIAL = {
        0: "nt",
        1 * EMUL: "g",
        (F_LEPTON * FMUL) + (+1 * ZMUL): "e+",
        (F_LEPTON * FMUL) + (-1 * ZMUL): "e-",
        (F_LEPTON * FMUL) + (+1 * EMUL): "nue",
        (F_LEPTON * FMUL) + (+2 * EMUL): "nueb",
        VOID_IDX: VOID_STRING,
    }

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
        void=True,
    ):
        """
        Initialize Isotope.

        isomer == True:
            interpret istopes as gs isomers
        """
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
            elif B is None:
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
                self.E == 0
        s = name
        if self._is_ion(s):
            self.B, self.F, self.Z, self.A, self.E = s.tuple()
            # conversion to Ion class
            # or could add keyword
            if not issubclass(type(self), type(s)):
                self.B = 0
                self.F = 0
            s = None

        if idx is None:
            try:
                i = int(s)
            except (ValueError, TypeError):
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
            if len(s) == 1:
                self.Z = s[0]
                self.F = self.F_ELEMENT
            if len(s) == 2:
                self.Z = s[0]
                self.A = s[1]
                self.F = self.F_ISOTOPE
            if len(s) == 3:
                self.Z = s[0]
                self.A = s[1]
                self.E = s[2]
                self.F = self.F_ISOMER
            if len(s) == 4:
                self.F = s[0]
                self.Z = s[1]
                self.A = s[2]
                self.E = s[3]
            if len(s) == 5:
                self.B = s[0]
                self.F = s[1]
                self.Z = s[2]
                self.A = s[3]
                self.E = s[4]
            s = None
        if isinstance(s, Ion):
            self.B, self.F, self.Z, self.A, self.E = s.tuple()
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
        if self.F in (self.F_ISOTOPE, self.F_ISOMER):
            if self.A < self.Z:
                self.idx = self.VOID_IDX
                self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(self.idx)
        if self.idx == self.VOID_IDX and not void:
            raise AttributeError("Could not create valid isotope from provided input.")
        self.N = self.get_N()
        if self.__class__ == Ion:
            super().__init__()

    @classmethod
    def _is_ion(cls, ion):
        return isinstance(ion, cls) or cls.__name__ == type(ion).__name__

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

        if z > 0 and idx < 0:
            z = -z
        elif idx == cls.VOID_IDX:
            b, f, z, a, e = cls.VOID
        elif idx < 0:
            e = -e
        return b, f, z, a, e

    @classmethod
    def ion2idx(cls, s=""):
        return cls.bfzae2idx(cls.ion2bfzae(s))

    def update_idx(self):
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

    def get_N(self):
        N = self.A - self.Z
        if self.F in (
            self.F_BOSON,
            self.F_ELEMENT,
            self.F_ISOBAR,
            self.F_HADRON,
            self.F_LEPTON,
        ):
            N = 0
        return N

    def tuple(self):
        return (self.B, self.F, self.Z, self.A, self.E)

    def dict(self):
        return {"F": self.F, "Z": self.Z, "A": self.A, "E": self.E, "B": self.B}

    def index(self):
        return self.idx

    def index_ZAE(self):
        return self.E * self.EMUL + self.Z * self.ZMUL + self.A * self.AMUL

    def ZAE(self):
        return (
            self.Z,
            self.A,
            self.E,
        )

    # A,Z, ... should become properties
    # def get_idx(self):
    #     return self.idx

    # def get_Z(self):
    #     return self.Z

    # def get_A(self):
    #     return self.A

    # def get_E(self):
    #     return self.E

    # def get_F(self):
    #     return self.F

    # def get_B(self):
    #     return self.B

    def Ye(self):
        if self.Z < 0 or self.A == 0:
            return -1
        return self.Z / self.A

    def mu(self):
        if self.Z < 0:
            return -1
        return self.A / (self.Z + 1)

    def mue(self):
        if self.Z <= 0:
            return -1
        return self.A / self.Z

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

    # def element(self):
    #     return Ion(Z = self.Z)

    # def isobar(self):
    #     return Ion(A = self.A)

    # def isotone(self):
    #     return Ion(N = self.N)

    # def isotope(self):
    #     return Ion(A = self.A,
    #                Z = self.Z)

    def element_symbol(self, upcase=True):
        s = ""
        if self.Z >= 0:
            if upcase:
                s = Elements[self.Z]
            else:
                s = elements[self.Z]
        return s

    def LaTeX(self, math=False):
        if self.is_element():
            el = Elements[self.Z]
            if math:
                return r"\mathrm{" + el + "}"
            else:
                return el
        a = str(self.A)
        el = Elements[self.Z]
        if math:
            return r"^{" + a + r"}\mathrm{" + el + "}"
        else:
            return r"$^{" + a + "}$" + el

    def NuGrid(self):
        el = Elements[self.Z]
        if self.is_element():
            return el
        A = str(self.A)
        if self.is_isotope():
            return f"{el}-{A}"
        if self.is_isomer():
            return f"{el}-{A}{self.isomer_name(self.E, m1=True)}"
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
        else:
            return f"{cls.ISOMER:s}{E:d}"

    def mpl(self):
        if self.is_element():
            el = Elements[self.Z]
            return r"$\mathsf{{{:}}}$".format(el)
        if self.is_isobar():
            return r"$\mathsf{{{:d}}}$".format(self.A)
        if self.is_isotone():
            return r"$\mathsf{{{:d}}}$".format(self.N)
        if self.is_isotope():
            el = Elements[self.Z]
            a = str(self.A)
            return r"$\mathsf{{^{{{:}\!}}{:}}}$".format(a, el)
        if self.is_isomer():
            el = Elements[self.Z]
            a = str(self.A)
            m = self.isomer_name(self.E)
            return r"$\mathsf{{^{{{:}\!}}{:}^{{{:}\!}}}}$".format(a, el, m)
        return "*"

    def LaTeX_Table(
        self, math=False, nozero=False, protons=False, manual=False, width=11, align=0
    ):
        a = str(self.A)
        el = Element[self.Z]
        if manual:
            s = "$^{" + a + "}$&" + el
        else:
            s = r"\I{" + a + "}{" + el + "}"
        if nozero and self.idx == self.VOID_IDX:
            s = "--"
        if align == -1:
            s = s.ljust(width)
        elif align == 1:
            s = s.rjust(width)
        else:
            s = s.center(width)
        if math and not manual:
            s = "$" + s + "$"
        return s

    def _name(self, upcase=True):
        """
        Convert normal isotopes and isomers to string.

        For isomers use 'g', 'm', 'm2', m3', ...
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
        s = self.SPECIAL.get(self.idx, s)
        return s

    @classmethod
    def ion2bfzae(cls, s="", element=False, isomer=False):
        sz, sa, se = cls._decompose(s, element=element)
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
        if sz in ("z=", "Z=", "Z:", "Z:", "z-", "Z-", "z", "Z"):
            valid = a > 0
            z = a
            a = 0
            try:
                z = int(sa)
            except ValueError:
                valid = False
            f = cls.F_ELEMENT
            se = 0
        try:
            z = Elements.index(sz)
        except ValueError:
            pass
        if z is None:
            try:
                z = elements.index(sz)
            except ValueError:
                pass
        if z is None:
            try:
                z = z2Name.index(sz.strip("-"))
            except ValueError:
                pass
        if z is None:
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
        elif sz == "e-" and sa == "" and se == "":
            z = -1
            e = 0
            f = cls.F_LEPTON
            valid = True
        elif (sz == "e+") and sa == "" and se == "":
            z = +1
            e = 0
            f = cls.F_LEPTON
            valid = True
        elif s == "g":
            # todo - add spin/E level
            f = cls.F_BOSON
            valid = True
        elif s == "nue":
            f = cls.F_LEPTON
            e = 1
            z = 0
            valid = True
        elif s == "nueb":
            f = cls.F_LEPTON
            e = 2
            z = 0
            valid = True
        if f == cls.F_ISOTOPE and z > 0 and a == 0:
            f = cls.F_ELEMENT
        if f == cls.F_ISOMER and z > 0 and a == 0:
            valid = False
        if not valid:
            b, f, z, a, e = cls.VOID
        return b, f, z, a, e

    def isomer(self, Z=None, A=None, N=None, E=map):
        """
        return self if isomer, isomer if enough info is provided, or VOID otherwise

        """
        if self.is_isomer() and E is None:
            if A is None and N is None and Z is None:
                return self
            E = self.E
        ion = self.isotope(A=A, Z=Z, N=N)
        if ion is VOID:
            return ion
        if E == "map" or E is map:
            E = isomermap
        if isinstance(E, Mapping):
            try:
                E = E[ion]
            except KeyError:
                return VOID
        if E is None:
            return VOID
        return self.__class__(A=ion.A, Z=ion.Z, E=E)

    def isotope(self, Z=None, A=None, N=None):
        """
        return isotope if isomer, or self if isotope, or VOID otherwise
        """
        if self.is_isotope() and Z == A == N is None:
            return self

        n = 3 - [Z, A, N].count(None)
        if self.is_isotope() or self.is_isomer():
            if n < 2 and Z is None:
                Z = self.Z
                n += 1
            if n < 2 and A is None:
                A = self.A
                n += 1
        elif self.is_isobar():  # has A
            if n < 2 and A is None:
                A = self.A
        elif self.is_isotone():  # has N
            if n < 2 and N is None:
                N = self.N
        elif self.is_element():  # has Z
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
                    return Ion(N=N, A=A)
        else:
            if A is None:
                if N is None:
                    return VOID
                else:
                    return Ion(Z=Z, N=N)
            else:
                if N is None:
                    if Z > A:
                        return VOID
                    return Ion(Z=Z, A=A)
                else:
                    if N + Z == A:
                        return Ion(Z=Z, A=A)
        return VOID

    def isobar(self, A=None):
        """
        return isobar if isomer, isotope, or self if isobar, or VOID otherwise
        """
        if A is not None:
            return Ion(A=A)
        if self.is_isobar():
            return self
        if self.is_isotope() or self.is_isomer:
            return Ion(A=self.A)
        return self.VOID

    def isotone(self, N=None):
        """
        return isotone if isomer, or isotope, or self if isotone, or VOID otherwise
        """
        if N is not None:
            return Ion(N=N)
        if self.is_isotone():
            return self
        if self.is_isotope() or self.is_isomer:
            return Ion(N=self.N)
        return self.VOID

    def element(self, Z=None):
        """
        return element if isomer, or isotope, or self if element, or VOID otherwise
        """
        if Z is not None:
            return Ion(Z=Z)
        if self.is_element():
            return self
        if self.is_isotope() or self.is_isomer:
            return Ion(Z=self.Z)
        return self.VOID

    @classmethod
    def _decompose(cls, s="", element=False):
        """
        This routine TRIES to decompose string in
        1) element symbol
        2) mass number
        3) excited state number
        All are returned as strings.

        The optional parameter 'element' enforces that
        'p' is taken as *element* potassium (not H1)
        'n' is taken as *element* nitrogen (not neutron)
        """
        name = s.strip()
        n = len(name)
        el = ""
        a = ""
        e = ""

        # get numbers
        n = re.findall(r"\d+", name)

        # get strings
        cx = re.findall(r"\D+", name)

        c = []
        for x in cx:
            xx = x.split("-")
            cy = [y for y in xx if y != ""]
            c += cy
        if len(c) == 2:
            if c[0] in ("m", "g"):
                c = c[::-1]
            if c[0][0] == "*":
                c = c[::-1]
        if len(n) > 0:
            a = n[0]
        if len(n) > 1:
            e = n[1]
        if len(n) > 2:
            raise ValueError(f"Can't understand isotope '{s}'.")
        if len(c) > 0:
            el = c[0]
        if len(el) > 0:
            if el[-1] in cls.EXCITE and len(c) == 1 and len(n) == 2:
                c.append(el[-1])
                el = el[:-1]
        if len(c) == 2:
            if c[1] == "g":
                e = "0"
            if c[1] == "m" and len(n) == 1:
                e = "1"
            if c[1][0] == "*" and len(n) == 1:
                e = str(len(c[1]))
                assert c[1].count("*") == len(c[1])
            if not c[1] in ("m", "g") and not c[1][0] == "*":
                raise ValueError(f"Can't understand isotope '{s}'.")

        if len(c) == 1 and c[0][-1] == "*":
            e = 0
            while c[0][-1] == "*":
                c[0] = c[0][:-1]
                e += 1
            e = str(e)
            el = c[0]

        if len(c) == 1 and c[0][0] == "*":
            e = 0
            while c[0][0] == "*":
                c[0] = c[0][1:]
                e += 1
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
            a = ""
        elif el in ("d", "D"):
            el = "H"
            if a != "":
                raise AttributeError('"d" already implies mass')
            a = "2"
        elif el in ("t", "T"):
            el = "H"
            if a != "":
                raise AttributeError('"t" already implies mass')
            a = "3"
        elif (element) and s == "n":
            el = "N"
        elif s == "n":
            el = "nt"
            a == 1
        elif el in ("n", "nt") and a == "1":
            el = "nt"
        elif s == "g":
            el = ""
            a = ""
            e = "1"
        elif s.lower() in ("e-", "b-", "bd", "pc"):
            el = "e-"
        elif (s.lower() in ("e+", "b+", "ec")) or (
            (not elements) and (s.lower() == "pd")
        ):
            el = "e+"
        el = el.strip()
        a = a.strip()
        e = e.strip()
        return el, a, e

    def dexcite(self):
        self.E = 0
        self.idx = self.bfzae2idx(F=self.F, Z=self.Z, A=self.A, E=self.E, B=self.B)

    def dexcited(self):
        return self.__class__(
            F=self.F,
            A=self.A,
            Z=self.Z,
            B=self.B,
        )

    def __str__(self):
        return self._name(upcase=True)

    def _add(self, x, sign1=+1, sign2=+1):
        if not self._is_ion(x):
            y = self.__class__(x)
            if y.idx == self.VOID_IDX:
                raise AttributeError(f"Cannot convert {x!r} to {self.__class__!r}.")
            else:
                x = y
        A = sign1 * self.A + sign2 * x.A
        Z = sign1 * self.Z + sign2 * x.Z
        E = sign1 * self.E + sign2 * x.E
        if not (
            ((self.F & self.F_GROUP_MASK == self.F_ISOMER) and x.is_photon())
            or (x.F & self.F_GROUP_MASK == self.F_ISOMER)
            and x.is_photon()
        ):
            E = None
        return self.__class__(
            Z=Z,
            A=A,
            E=E,
        )

    def __add__(self, x):
        return self._add(x)

    __radd__ = __add__

    def __sub__(self, x):
        return self._add(x, +1, -1)

    def __rsub__(self, x):
        return self._add(x, -1, +1)

    def __mul__(self, x):
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
        if not isinstance(x, Ion):
            x = Ion(x)
        return self.idx < x.idx

    def __eq__(self, x):
        if not isinstance(x, Ion):
            x = Ion(x)
        return self.idx == x.idx

    def __len__(self):
        return 1

    def __repr__(self):
        base = f"('{self.Name():s}')"
        # if self.__class__ == Ion:
        #     if self.is_isomer():
        #         return 'Isomer' + base
        #     if self.is_isotope():
        #         return 'Isotope' + base
        #     if self.is_isobar():
        #         return 'Isobar' + base
        #     if self.is_isotine():
        #         return 'Isotone' + base
        #     if self.is_element():
        #         return 'Element' + base
        #     if self.is_lepton():
        #         return 'Lepton' + base
        #     if self.is_hadron():
        #         return 'Hadron' + base
        return self.__class__.__name__ + base

    def __getstate__(self):
        # need to add version number
        # should save index instead
        #        print('xxx_save')
        return self.tuple()

    def __setstate__(self, x):
        # need to add check of version number
        #        print('xxx_load')
        self.__init__(x)

    # this does not appear to be called for a new-style class ever.
    # def __getinitargs__(self):
    #     print 'xxx_old'
    #     return (self.tuple(),)
    # the following could work as well, but does appear less flexible
    # def __getnewargs__(self):
    #     print 'xxx_new'
    #     return (self.tuple(),)
    def element_name(self):
        return z2name[self.Z]

    def is_nucleus(self):
        return self.is_isomer() or self.is_isotope()

    def is_lepton(self):
        return self.F & self.F_GROUP_MASK == self.F_LEPTON

    def is_isomer(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOMER

    def is_isotope(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTOPE

    def is_isobar(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOBAR

    def is_isotone(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTONE

    def is_element(self):
        return self.F & self.F_GROUP_MASK == self.F_ELEMENT

    def is_hadron(self):
        return self.F & self.F_GROUP_MASK == self.F_HADRON

    def is_boson(self):
        return self.F & self.F_GROUP_MASK == self.F_BOSON

    def is_photon(self):
        return (
            (self.F & self.F_GROUP_MASK == self.F_BOSON)
            and (self.Z == 0)
            and (self.A == 0)
            and (self.E == 1)
        )

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
        np.frompyfunc(lambda x: x.is_lepton(), 1, 1)(y), dtype=bool
    )
    ufunc_is_isotope = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_isotope(), 1, 1)(y), dtype=bool
    )
    ufunc_is_isobar = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_isobar(), 1, 1)(y), dtype=bool
    )
    ufunc_is_isotone = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_isotone(), 1, 1)(y), dtype=bool
    )
    ufunc_is_element = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_element(), 1, 1)(y), dtype=bool
    )
    ufunc_is_isomer = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_isomer(), 1, 1)(y), dtype=bool
    )
    ufunc_is_hadron = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_hadron(), 1, 1)(y), dtype=bool
    )
    ufunc_is_nucleus = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_nucleus(), 1, 1)(y), dtype=bool
    )
    ufunc_is_photon = lambda y: np.array(
        np.frompyfunc(lambda x: x.is_photon(), 1, 1)(y), dtype=bool
    )

    ufunc_type = lambda y: np.array(
        np.frompyfunc(lambda x: x.type, 1, 1)(y), dtype=np.str
    )

    ufunc_idx = lambda y: np.array(np.frompyfunc(lambda x: x.idx, 1, 1)(y), dtype=int)

    ufunc_ion = lambda y: np.array(
        np.frompyfunc(lambda x: Ion(x), 1, 1)(y), dtype=object
    )
    ufunc_ion_from_idx = lambda y: np.array(
        np.frompyfunc(lambda x: Ion(idx=x), 1, 1)(y), dtype=object
    )

    @CachedAttribute
    def type(self):
        if self.is_lepton():
            return "lepton"
        if self.is_isotope():
            return "isotope"
        if self.is_isobar():
            return "isobar"
        if self.is_isotone():
            return "isotone"
        if self.is_element():
            return "element"
        if self.is_isomer():
            return "isomer"
        if self.is_hadron():
            return "hadron"
        if self.is_nucleon():
            return "nucleon"
        if self.is_photon():
            return "photon"
        if self.is_boson():
            return "boson"
        return "unkown"

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
                return self()
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
        Ion('c12')('p','g') --> Ion('N13')

        Ion('c12')('p') -->

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


GAMMA = Ion(Ion.GAMMA)
VOID = Ion(Ion.VOID)
NEUTRON = Ion(Z=0, A=1)
PROTON = Ion(Z=1, A=1)
ALPHA = Ion(Z=2, A=4)


# some convenience functions
def ionarr(arr):
    if isinstance(arr, str):
        arr = re.split("[^-0-9a-zA-Z*]+", arr)
    if not isinstance(arr, Iterable):
        arr = (arr,)
    return np.array([Ion(a) for a in arr], dtype=object)


# the following classes need to be fleshed out and used consequnetly
class Element(Ion):
    pass


class Isotone(Ion):
    pass


class Isobar(Ion):
    pass


class Isotope(Ion):
    pass


class Isomer(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs["isomer"] = True
        super().__init__(*args, **kwargs)


class Photon(Ion):
    pass


class Lepton(Ion):
    pass


class Hadron(Ion):
    pass


# create functions
def element(Z):
    return Ion(Z=Z)


def isotone(N):
    return Ion(N=N)


def isobar(A):
    return Ion(A=A)


def isotope(A, Z):
    return Ion(A=A, Z=Z)


def isomer(A, Z, E=0):
    return Ion(A=A, Z=Z, E=E)


def photon(E=1):
    return Ion(E=E)


def lepton(Z=0, E=1):
    return Ion(Z=Z, E=E, lepton=True)


def hadron(Z=0, E=1):
    return Ion(Z=Z, E=E, F=Ion.F_HARDON)


class IonFactory:
    """
    Factory function class.

    This one should be used to allow for later refactoring in derived
    classes for different ion types.

    TODO - implement factory function

    """

    def __call__(self, *args, **kwargs):
        ion = Ion(*args, **kwargs)
        if ion.is_isomer():
            return Isomer(ion)
        if ion.is_isotope():
            return Isotope(ion)
        if ion.is_element():
            return Element(ion)
        if ion.is_isobar():
            return Isobar(ion)
        if ion.is_isotone():
            return Isotone(ion)
        if ion.is_photon():
            return Photon(ion)
        if ion.is_lepton():
            return Lepton(ion)
        if ion.is_hadron():
            return Hadron(ion)
        return ion

    def __getattr__(self, attr):
        x = Ion(attr)
        if x is not VOID:
            return x
        raise AttributeError(attr)


ion = IonFactory()


class IsomerMap(collections.defaultdict):
    """
    define map that can be used as 'E' paramater in ion.isomer function

    TODO - add more isotopes, provide function to load from file
    """

    def __init__(self, default=lambda: 0, map=None):
        super().__init__(default)
        if map is not None:
            self.update(map)

    def __getitem__(self, item):
        if not isinstance(item, Ion):
            item = Ion(item)
        return super().__getitem__(item)

    def __setitem__(self, item, value):
        if not isinstance(item, Ion):
            item = Ion(item)
        super().__setitem__(item, value)


isomermap = IsomerMap(map={Ion("ta180"): 1})

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
        assert bit not in other_bits_register, "bit already exisits"

    other_bits_register[bit] = (name, class_)
