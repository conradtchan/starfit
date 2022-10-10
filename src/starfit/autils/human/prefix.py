"""
Defines prefix class and default prefixes

Unofficial:
   zepto (1e-21)
   yocto (1e-24)
"""

from copy import copy

from numpy import isreal


class Prefix(object):
    def __init__(
        self,
        name,
        base,
        power=1,
        *,
        latex=None,
        unicode=None,
        output=True,
    ):
        self._name = name
        self._base = base
        self._power = power
        self._latex = latex
        self._unicode = unicode
        self._output = output

    @property
    def output(self):
        return self._output

    @property
    def value(self):
        if self._power == 0:
            return 1
        elif self._power < 0:
            return 1 / self._base**-self._power
        return self._base**self._power

    @property
    def name(self):
        if self._name is None:
            return self._unicode
        return self._name

    @property
    def latex(self):
        if self._latex is None:
            return self.name
        return self._latex

    @property
    def unicode(self):
        if self._unicode is None:
            return self.name
        return self._unicode

    def __pow__(self, other):
        if isreal(other):
            new = copy(self)
            new._power += other
            return new
        return NotImplemented()

    def __mul__(self, other):
        if hasattr(other, "value"):
            other = getattr(other, "value")
        if isreal(other):
            if self._power == 0:
                return other
            elif self._power < 0:
                return other / self._base**-self._power
            return other * self._base**self._power
        return NotImplemented()

    def __rtruediv__(self, other):
        if hasattr(other, "value"):
            other = getattr(other, "value")
        if isreal(other):
            if self._power == 0:
                return other
            elif self._power < 0:
                return other * self._base**-self._power
            return other / self._base**self._power
        return NotImplemented()

    def __repr__(self):
        if self._name == "":
            name = ""
        else:
            name = f"{self.name}, "
        if self._power == 0:
            base = "1"
        else:
            base = f"{self._base}"
        if self._power in (0, 1):
            power = ""
        else:
            power = f"**{self._power}"
        if self._power == 0 and name == "":
            base = ""
        return f"{self.__class__.__name__}({name}{base}{power})"


Prefix.neutral = Prefix("", 1, 0)

_Prefixes = (
    Prefix("k", 10, 3),
    Prefix("M", 10, 6),
    Prefix("G", 10, 9),
    Prefix("T", 10, 12),
    Prefix("P", 10, 15),
    Prefix("E", 10, 18),
    Prefix("Z", 10, 21),
    Prefix("Y", 10, 24),
)

_Prefixes2 = (
    Prefix("da", 10, 1, output=False),
    Prefix("h", 10, 2, output=False),
    Prefix("k", 10, 3),
    Prefix("M", 10, 6),
    Prefix("G", 10, 9),
    Prefix("T", 10, 12),
    Prefix("P", 10, 15),
    Prefix("E", 10, 18),
    Prefix("Z", 10, 21),
    Prefix("Y", 10, 24),
)

_Prefixes3 = (
    Prefix("da", 10, 1),
    Prefix("h", 10, 2),
    Prefix("k", 10, 3),
    Prefix("M", 10, 6),
    Prefix("G", 10, 9),
    Prefix("T", 10, 12),
    Prefix("P", 10, 15),
    Prefix("E", 10, 18),
    Prefix("Z", 10, 21),
    Prefix("Y", 10, 24),
)

_prefixes = (
    Prefix("m", 10, -3),
    Prefix("u", 10, -6, unicode="\u03BC", latex=r"{\mu}"),
    Prefix("n", 10, -9),
    Prefix("p", 10, -12),
    Prefix("f", 10, -15),
    Prefix("a", 10, -18),
    Prefix("z", 10, -21),
    Prefix("y", 10, -24),
)

_prefixes2 = (
    Prefix("d", 10, -1, output=False),
    Prefix("c", 10, -2, output=False),
    Prefix("m", 10, -3),
    Prefix("u", 10, -6, unicode="\u03BC", latex=r"{\mu}"),
    Prefix("n", 10, -9),
    Prefix("p", 10, -12),
    Prefix("f", 10, -15),
    Prefix("a", 10, -18),
    Prefix("z", 10, -21),
    Prefix("y", 10, -24),
)

_prefixes3 = (
    Prefix("d", 10, -1),
    Prefix("c", 10, -2),
    Prefix("m", 10, -3),
    Prefix("u", 10, -6, unicode="\u03BC", latex=r"{\mu}"),
    Prefix("n", 10, -9),
    Prefix("p", 10, -12),
    Prefix("f", 10, -15),
    Prefix("a", 10, -18),
    Prefix("z", 10, -21),
    Prefix("y", 10, -24),
)

_BinPrefixes = (
    Prefix("ki", 2, 10),
    Prefix("Mi", 2, 20),
    Prefix("Gi", 2, 30),
    Prefix("Ti", 2, 40),
    Prefix("Pi", 2, 50),
    Prefix("Ei", 2, 60),
    Prefix("Zi", 2, 70),
    Prefix("Yi", 2, 80),
)

_binPrefixes = (
    Prefix("mi", 2, -10),
    Prefix("ui", 2, -20, unicode="\u03BCi", latex=r"{\mu}i"),
    Prefix("ni", 2, -30),
    Prefix("pi", 2, -40),
    Prefix("fi", 2, -50),
    Prefix("ai", 2, -60),
    Prefix("zi", 2, -70),
    Prefix("yi", 2, -80),
)
