from copy import copy
from math import inf, nan

from .power import Power
from .prefix import Prefix
from .unit import Unit


# keep parts in braces in LaTeX math mode and rest is put into mathrm mode
def convert_latex(su):
    sx = ""
    lc = 0
    cc = 0
    for sy in su:
        if sy == "{":
            if lc > 0:
                sx += sy
            elif cc > 0:
                sx += "}"
                cc = 0
            lc += 1
        elif sy == "}":
            if lc > 1:
                sx += sy
            lc -= 1
        elif lc == 0:
            if cc == 0:
                cc += 1
                sx += r"\mathrm{"
            sx += sy
        else:
            sx += sy
    else:
        assert lc == 0
        if cc > 0:
            sx += "}"
            cc = 0
    return sx


_latex_sep = r"{\,}"
_text_sep = r" "
_unicode_sep = "\u2009"


class Range(object):
    def __init__(self, *args, **kwargs):
        args_ = list()
        start = kwargs.pop("start", None)
        stop = kwargs.pop("stop", None)
        unit = kwargs.pop("unit", None)
        prefix = kwargs.pop("prefix", None)
        power = kwargs.pop("power", None)
        exponent = kwargs.pop("exponent", None)
        if len(kwargs) > 0:
            raise AttributeError(f"unexpected keywords: {kwargs}")

        for arg in args:
            if isinstance(arg, Unit):
                assert unit is None
                unit = arg
                continue
            if isinstance(arg, Prefix):
                assert prefix is None
                prefix = arg
                continue
            if isinstance(arg, Power):
                assert power is None
                power = arg
                continue
            if isinstance(arg, bool):
                assert exponent is None
                exponent = arg
                continue
            args_.append(arg)
        if unit is None:
            unit = Unit.neutral
        if prefix is None:
            prefix = Prefix.neutral
        if power in (
            None,
            Ellipsis,
            nan,
        ):
            power = Power.neutral
        elif not isinstance(power, Power):
            # raise NotImplementedError(f'power = {power}')
            power = Power(value=Power)
        if exponent is None:
            exponent = False
        if len(args_) > 0 and start is None:
            start = args_.pop(0)
        if start is None:
            start = -inf
        if len(args_) > 0 and stop is None:
            stop = args_.pop(0)
        if len(args_) > 0:
            raise AttributeError(f"unexpected arguments: {args_}")

        if stop is None:
            stop = +inf
        if start > stop:
            stop, start = start, stop
        self._start = start
        self._stop = stop
        self._unit = unit
        self._power = power
        self._prefix = prefix
        self._unit = unit
        self._exponent = exponent

    @property
    def value(self):
        return self._unit.value * self._prefix.value * self._power.value

    def str(self, unicode=False, latex=False, unit=False):
        if unit:
            if unicode:
                return self.unicode_unit
            if latex:
                return self.latex_unit
            return self.name_unit
        if unicode:
            return self.unicode
        if latex:
            return self.latex
        return self.name

    @property
    def name(self):
        e = self._power.exponent()
        u = self._prefix.name + self._unit.name
        if len(u) >= 0:
            return e + _text_sep + u
        return e

    @property
    def name_unit(self):
        e = self._power.unit(fill=_text_sep)
        return e + self._prefix.name + self._unit.name

    @property
    def unicode(self):
        e = self._power.exponent_unitcode(fill=_unicode_sep)
        return e + _unicode_sep + self._prefix.unicode + self._unit.unicode

    @property
    def unicode_unit(self):
        e = self._power.unit_unitcode(fill=_unicode_sep)
        return e + _unicode_sep + self._prefix.unicode + self._unit.unicode

    @property
    def latex(self):
        e = self._power.exponent_latex(fill=_latex_sep)
        su = e + _latex_sep + self._prefix.latex + self._unit.latex
        return convert_latex(su)

    @property
    def latex_unit(self):
        e = self._power.unit_latex(fill=_latex_sep)
        su = e + _latex_sep + self._prefix.latex + self._unit.latex
        return convert_latex(su)

    @property
    def prefix_name(self):
        return self._unit.prefix

    @property
    def prefix_value(self):
        return self._prefix.value

    @property
    def unit_name(self):
        return self._unit.name

    @property
    def unit_value(self):
        return self._unit.value

    @property
    def power_name(self):
        return self._power.name

    @property
    def power_value(self):
        return self._power.value

    @property
    def prefix(self):
        return self._prefix

    @property
    def has_prefix(self):
        return (self._prefix.name != "") or (self._prefix.value != 1)

    @property
    def unit(self):
        return self._unit

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def exponent(self):
        return self._exponent

    @property
    def power(self):
        return self._power

    def set_start(self, start):
        self._start = start

    def add_power(self, power):
        power = self._power + power
        new = copy(self)
        new._power = power
        return new

    @property
    def pure(self):
        if self.exponent is True:
            return False
        if self.power.value != 1:
            return False
        if self.prefix.value != 1:
            return False
        return True

    def __add__(self, other):
        if isinstance(other, Power):
            return self.add_power(other)
        raise AttributeError()

    def __repr__(self):
        return f"{self.__class__.__name__}({self._start}, {self._stop}, {self._unit}, {self._prefix}, {self._power!r}, {self._exponent})"


def Ranges(object):
    def __init__(self):
        self._data = list()

    def __getitem__(self, key):
        return self.__data[key]

    def insert(self, item):
        pass
