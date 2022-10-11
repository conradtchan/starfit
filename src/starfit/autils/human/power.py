from math import log10


class Power(object):
    def __init__(self, power=None, value=None):
        if power is None and value is None:
            self._power = 0
            self._value = 1
        elif value is not None and power is not None:
            raise AttributeError("can only set power or value")
        elif value is not None:
            if isinstance(value, str):
                value = (
                    value.strip().replace("d", "e").replace("D", "e").replace("E", "e")
                )
                if value.startswith("e"):
                    value = "1" + value
                if value == "":
                    value = 1
            self._value = float(value)
            self._power = int(round(log10(self._value)))
            assert abs(10**self._power / self._value - 1) < 1e-12
        else:
            self._power = int(power)
            self._value = 10**self._power

    def __str__(self):
        if self._power == 0:
            return "1"
        return f"1e{self._power:d}"

    @property
    def value(self):
        return self._value

    @property
    def power(self):
        return self._power

    def unit(self, latex=False, unicode=False, fill=""):
        if latex:
            return self.unit_latex(fill=fill)
        if unicode:
            return self.unit_unicode(fill=fill)
        if self._power == 0:
            return ""
        return f"1e{self._power:d}" + fill

    def exponent(self, latex=False, unicode=False, fill="", unit=False):
        if unit:
            return self.unit(latex, unicode, fill)
        if latex:
            return self.exponent_latex(fill)
        if unicode:
            return self.exponent_unicode()
        if self._power == 0:
            return ""
        return f"e{self._power:d}" + fill

    def exponent_latex(self, fill=""):
        if self._power == 0:
            return ""
        if self._power == 1:
            return "{\times10}" + fill
        return rf"{{\times10^{{{self._power:d}}}}}" + fill

    def unit_latex(self, fill=""):
        if self._power == 0:
            return ""
        if self._power == 1:
            return "{10}" + fill
        return rf"{{10^{{{self._power:d}}}}}" + fill

    _power_unicode_map = {
        "0": "\u2070",
        "1": "\u00b9",
        "2": "\u00b2",
        "3": "\u00b3",
        "4": "\u2074",
        "5": "\u2075",
        "6": "\u2076",
        "7": "\u2077",
        "8": "\u2078",
        "9": "\u2079",
        "+": "\u207a",
        "-": "\u207b",
    }

    def exponent_unicode(self, fill=""):
        if self._power == 0:
            return ""
        if self._power == 1:
            return "\u22c510" + fill
        s = f"{self._power:d}"
        pow = str.join("", (self._power_unicode_map[x] for x in s))
        return "\u22c510" + pow + fill

    def unit_unicode(self, fill=""):
        if self._power == 0:
            return ""
        if self._power == 1:
            return "10" + fill
        s = f"{self._power:d}"
        pow = str.join("", (self._power_unicode_map[x] for x in s))
        return "10" + pow + fill

    def __repr__(self):
        return f"{self.__class__.__name__}({self._power:d})"

    def __add__(self, other):
        if isinstance(other, self.__class__):
            power = self._power + other._power
            return self.__class__(power)
        raise AttributeError()


Power.neutral = Power()
