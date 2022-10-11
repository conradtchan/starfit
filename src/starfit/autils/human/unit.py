import numpy as np


class Unit(object):
    def __init__(
        self,
        name=None,
        value=None,
        *,
        div_lim=None,
        latex=None,
        unicode=None,
        prefix=...,
        output=True,
        out_prefix=True,
        power_limit=None,  # values relative to unit value
        hi_power_limit=...,  # ... default
        lo_power_limit=...,
        prefix_power=True,  # True | False
        lo_prefix_power=None,
        hi_prefix_power=None,
        prefix_dimension=1,
        alt=None,
    ):
        if name in (
            None,
            Ellipsis,
        ):
            name = ""
        if not isinstance(name, str):
            name = list(name)
        if isinstance(name, list):
            name, *alt_ = name
        else:
            alt_ = []
        self.name = name
        if alt is None:
            alt = []
        elif isinstance(alt, str):
            alt = [alt]
        else:
            alt = list(alt)
        self.alt = alt_ + alt
        if value in (
            None,
            Ellipsis,
        ):
            value = 1
        self.value = value
        self._latex = latex
        self._unicode = unicode
        self.div_lim = div_lim
        self.prefix = prefix
        self.output = output
        self.out_prefix = out_prefix
        self.prefix_dimension = prefix_dimension

        # limits for power
        if power_limit is not None:
            if isinstance(power_limit, str):
                power_limit = float(hi_power_limit)
            power_limit = np.exp(np.abs(np.log(power_limit)))
            if lo_power_limit is ...:
                lo_power_limit = 1 / power_limit
            if hi_power_limit is ...:
                hi_power_limit = power_limit
        else:
            if hi_power_limit is ...:
                hi_power_limit = None
            if lo_power_limit is ...:
                lo_power_limit = None
        if hi_power_limit is not None:
            if isinstance(hi_power_limit, str):
                hi_power_limit = float(hi_power_limit)
            if hi_power_limit > 0:
                hi_power_limit *= self.value
            else:
                hi_power_limit *= -1
        if lo_power_limit is not None:
            if isinstance(lo_power_limit, str):
                lo_power_limit = float(lo_power_limit)
            if lo_power_limit > 0:
                lo_power_limit *= self.value
            else:
                lo_power_limit *= -1
        self.hi_power_limit = hi_power_limit
        self.lo_power_limit = lo_power_limit

        # allow power on prefix
        if lo_prefix_power is None:
            lo_prefix_power = prefix_power
        if hi_prefix_power is None:
            hi_prefix_power = prefix_power
        self.lo_prefix_power = lo_prefix_power
        self.hi_prefix_power = hi_prefix_power

    @property
    def unicode(self):
        if self._unicode is None:
            return self.name
        return self._unicode

    @property
    def latex(self):
        if self._latex is None:
            return self.name
        return self._latex

    @property
    def names(self):
        return set([self.name, self.unicode]) | set(self.alt)

    def __repr__(self):
        if self.name == "":
            s = f"{self.value}"
        else:
            s = f"{self.name}, {self.value}"
        return f"{self.__class__.__name__}({s})"


Unit.neutral = Unit()


class Units(object):
    def __init__(self, *args):
        args = (Unit(**arg) if isinstance(arg, dict) else arg for arg in args)
        self.data = sorted(args, key=lambda x: x.value)
        self.lookup = {x.name: i for i, x in enumerate(self.data)}
        for i, x in enumerate(self.data):
            if x.value == 1:
                self.base = x
                self.base_index = i
                break
        else:
            raise Exception("[UNITS] No base unit found.")

    def __getitem__(self, index):
        if isinstance(index, str):
            return self.data[self.lookup[index]]
        return self.data[index]

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for x in self.data:
            yield x

    def output(self, base=None, large=None):
        # TODO - rewrite using index, to be set up in __init__
        for x in self.data:
            if x.output is False:
                continue
            if (base is not False or large is None) and x.value == 1:
                yield x
            elif large is not False and x.value > 1:
                yield x
            elif large is not True and x.value < 1:
                yield x

    def __repr__(self):
        s = ", ".join(repr(d) for d in self.data)
        return f"{self.__class__.__name__}({s})"
