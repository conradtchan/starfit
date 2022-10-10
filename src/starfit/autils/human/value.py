from math import ceil, floor, log10, nan

from numpy import isreal


class Value(object):
    def __init__(
        self,
        value,
        digits=None,
        decimals=None,
        comma=False,
        cut=True,
    ):
        if isinstance(value, str):
            try:
                value = int(value)
            except:
                pass
        if isinstance(value, str):
            try:
                value = float(value)
            except:
                pass
        ivalue = int(value)
        if value == ivalue:
            value = ivalue
            self._is_int = True
        else:
            self._is_int = False
        self._value = value
        self._comma = comma
        self._cut = cut
        self._digits = digits
        self._decimals = decimals

    @property
    def value(self):
        return self._value

    @property
    def float(self):
        return float(self)

    @property
    def name(self):
        return self.str()

    def str(self, digits=None, decimals=None, comma=None, cut=None):
        if digits is None:
            digits = self._digits
        if decimals is None:
            decimals = self._decimals
        if comma is None:
            comma = self._comma
        if cut is None:
            cut = self._cut
        if digits is not None:
            xval = self._value
            if abs(xval) > 1e-99:
                mag = int(ceil(log10(xval)))
            else:
                mag = 0
            decimals = digits - max(1, mag)  # leading 0 counts as digit
            if mag > digits:
                decimals += 1  # drop decimal point
            if round(xval, decimals) > 10**mag:
                decimals -= 1
            if decimals > 0:
                s = f"{xval:{digits+1}.{decimals}f}"
            else:
                s = f"{int(round(xval)):d}"
        elif decimals is None:
            s = f"{self.value:0f}"
            decimals = nan
        else:
            s = f"{self.value:0.{decimals}f}"
        if cut and s.count(".") > 0:
            s = s.rstrip("0").rstrip(".")
        if comma:
            l = len(s)
            j = s.find(".")
            if j == -1:
                j = l
            for i in range(j - 3, 0, -3):
                s = s[:i] + "," + s[i:]
        return s

    @property
    def is_int(self):
        return self._is_int

    @property
    def is_float(self):
        return not self._is_int

    @property
    def mag(self):
        return int(floor(log10(self.value)))

    @property
    def decimals(self):
        s = self.str()
        j = s.find(".")
        if j == -1:
            return 0
        return len(s) - j - 1

    def __mul__(self, other):
        if hasattr(other, "value"):
            other = getattr(other, "value")
        if isreal(other):
            return self.value * other
        return NotImplemented()

    def __len__(self):
        return len(str(self))

    def __str__(self):
        return self.str()

    def __repr__(self):
        return f"{self.__class__.__name__}({str(self)})"


Value.zero = Value(0)
Value.one = Value(1)
