debug = False

import re
from collections.abc import Iterable
from itertools import chain
from math import ceil, floor, inf, log10

from .config import Config
from .power import Power
from .prefix import Prefix
from .range import Range
from .unit import Unit
from .value import Value

epsf = 1 - 1 / (2 << 51)
epsm = 1 + 1 / (2 << 51)

# enums = (list, tuple, set, frozenset, )
enums = Iterable
negs = (
    "~",
    "!",
    "-",
)

hi_sentinels = {"+", ">"}
lo_sentinels = {"-", "<"}
any_sentinels = {..., True}
no_sentinels = {None, False}


class AmbiguousUnitError(Exception):
    def __init__(self, unit, prefix):
        name = prefix.name + unit.name
        super().__init__(
            f'Ambiguous or duplicate unit name or scale "{name}" '
            + f'from prefix "{prefix.name}" and base "{unit.name}".'
        )


# this may be good for float64 only
def _div_lim(x, digits=0):
    return x * (1 - 1 / (2 << 53)) - 0.5 * 10 ** (-digits)


class Human(object):
    """
    basis class for general scaling

    config has the following fields:

    convert:
       conversion function to basic unit type
    units:
       Units object

    TODO - it might be useful to also store upper limits at first
           and then consolidate later?
    """

    def __init__(self, config, /, **kwargs):
        self.debug = kwargs.get("debug", debug)

        if not isinstance(config, Config):
            config = Config(config)
        self.convert = config.convert
        self.integer = config.integer
        self.units = config.units
        self.hi_prefixes = config.hi_prefixes
        self.lo_prefixes = config.lo_prefixes
        self.zero_unit = config.zero_unit
        self.hi_prefix_base = config.hi_prefix_base
        self.lo_prefix_base = config.lo_prefix_base
        self.hi_power_base = config.hi_power_base
        self.lo_power_base = config.lo_power_base
        self.config = config

        # some defaults for later
        self.base = kwargs.get("base", 1)
        self.dec_lim = kwargs.get("dec_lim", 1)
        self.digits = kwargs.get("digits", 2)
        self.unicode = kwargs.get("unicode", False)

        # do some beautification of unit prefixes selection
        # replace indicators, convert all to sets
        #
        # exclude some units if requested on unit object as option
        # Ellipsis | True - allow all prefixes
        # None | False - no prefixes
        # '+', '>' - positive prefixes
        # '-', '<' - negative prefixes

        self.out_prefix = {""}
        self.out_prefix |= {p.name for p in self.lo_prefixes if p.output}
        self.out_prefix |= {p.name for p in self.hi_prefixes if p.output}
        self.out_prefix = frozenset(self.out_prefix)

        self.lo_prefixes_set = frozenset({p.name for p in self.lo_prefixes})
        self.hi_prefixes_set = frozenset({p.name for p in self.hi_prefixes})
        self.prefixes_set = frozenset(self.hi_prefixes_set | self.lo_prefixes_set)
        self.al_prefixes_set = frozenset(self.prefixes_set | set(("",)))
        for unit in self.units:
            for attr in (
                "prefix",
                "out_prefix",
            ):
                prefix = getattr(unit, attr)
                if isinstance(prefix, str):
                    prefix = {prefix}
                if isinstance(prefix, enums):
                    prefix = set(prefix)
                    if len(lo_sentinels & prefix) > 0:
                        prefix -= lo_sentinels
                        prefix |= self.lo_prefixes_set
                    if len(hi_sentinels & prefix) > 0:
                        prefix -= hi_sentinels
                        prefix |= self.hi_prefixes_set
                    if len(any_sentinels & prefix) > 0:
                        prefix -= any_sentinels
                        prefix |= self.prefixes_set
                    addall = None
                    for p in prefix:
                        if not isinstance(p, str):
                            addall = False
                            break
                        if p.startswith(negs):
                            addall = True
                        else:
                            addall = False
                            break
                    if addall is True:
                        prefix |= self.prefixes_set
                    # add all positives if there is only negations
                    addall = None
                    for p in prefix:
                        if not isinstance(p, str):
                            addall = False
                            break
                        if p.startswith(negs) and p[1:] in self.hi_prefixes_set:
                            addall = True
                        if p in self.hi_prefixes_set:
                            addall = False
                            break
                    if addall is True:
                        prefix |= self.hi_prefixes_set
                    # add all negatives if there is only negations
                    addall = None
                    for p in prefix:
                        if not isinstance(p, str):
                            addall = False
                            break
                        if p.startswith(negs) and p[1:] in self.lo_prefixes_set:
                            addall = True
                        if p in self.lo_prefixes_set:
                            addall = False
                            break
                    if addall is True:
                        prefix |= self.lo_prefixes_set
                    prefix_ = prefix.copy()
                    for p in prefix_:
                        if not isinstance(p, str):
                            continue
                        if p.startswith(negs):
                            if not p[1:] in prefix:
                                raise AttributeError(f'Invalid removal "{p}".')
                            prefix -= set(
                                (
                                    p,
                                    p[1:],
                                )
                            )
                elif prefix in lo_sentinels:
                    prefix = set(self.lo_prefixes_set)
                elif prefix in hi_sentinels:
                    prefix = set(self.hi_prefixes_set)
                elif prefix in any_sentinels:
                    prefix = set(self.prefixes_set)
                elif prefix in no_sentinels:
                    prefix = set()
                else:
                    raise AttributeError(f'Unknown selection for "{attr}": "{prefix}".')
                unknown = prefix - self.prefixes_set
                if len(unknown) > 0:
                    raise AttributeError(f'Unknown entries for "{attr}": "{unknown}".')
                setattr(unit, attr, prefix)
            unit.prefix = frozenset(unit.prefix)
            unit.out_prefix = frozenset(unit.prefix & unit.out_prefix & self.out_prefix)

            if self.debug:
                print(f"{unit.name}.unit.prefix: {unit.prefix}")
                print(f"{unit.name}.unit.out_prefix: {unit.out_prefix}")

        self._unit2scale = dict()
        for unit in self.units:
            for prefix in chain(self.lo_prefixes, (Prefix.neutral,), self.hi_prefixes):
                if (prefix.name not in unit.prefix) and (
                    prefix.name != Prefix.neutral.name
                ):
                    continue
                if unit.prefix_dimension != 1:
                    prefix **= unit.prefix_dimension
                scale = unit.value * prefix.value
                for unit_name in unit.names:
                    name = prefix.name + unit_name
                    if name in self._unit2scale:
                        raise AmbiguousUnitError(unit, prefix)
                    self._unit2scale[name] = (
                        scale,
                        unit,
                        prefix,
                    )

        self.outunits = list(self.units.output())

        # generate ranges for default settings
        self.ranges = dict()
        ranges_key = (self.digits, self.dec_lim)
        self.ranges[ranges_key] = self._compute_ranges(self.digits, self.dec_lim)

    def _compute_ranges(self, digits, dec_lim, /):

        divlim = 10 ** (1 + digits)

        def div_lim(scale):
            scale *= dec_lim
            mag = int(ceil(log10(scale)))
            decimals = digits - max(1, mag)  # leading 0 counts as digit
            if mag > digits:
                decimals += 1  # drop decimal point
            decimals += max(0, -int(floor(log10(dec_lim))))
            return _div_lim(scale, decimals)

        div_lim1 = div_lim(1)
        div_limd = div_lim(10**digits) / 10**digits

        # compute divlim scales.
        # set up output scalings as defined
        scales = []
        unit_ = None
        for unit in self.outunits:
            if unit_ is not None:
                if unit_.div_lim is not None:
                    scale = unit_.div_lim
                else:
                    scale = min(unit.value / unit_.value, divlim)
                    iscale = int(scale)
                    if iscale == scale:
                        scale = iscale
                scales.append((unit_, scale))
            unit_ = unit
        else:
            if unit.div_lim is not None:
                scale = unit.div_lim
            else:
                scale = divlim
            scales.append((unit, scale))

        # compute ranges
        ranges = list()

        unit, scale = scales[0]
        lim = unit.value * div_lim1  # preliminary, may be overwritten later
        prefix = Prefix.neutral
        exponent = False
        ranges.append(Range(lim, unit, prefix, exponent))
        unit_ = unit
        lim = unit.value * div_lim(scale)
        for iscale, (unit, scale) in enumerate(scales[1:]):
            lim0 = lim
            lim1 = unit.value * div_lim1
            # if unit.name == 'pc':
            #     breakpoint()
            if lim0 < lim1 * epsf:
                for junit, jscale in scales[iscale::-1]:
                    if not junit.output:
                        continue
                    prefixes = junit.out_prefix
                    for px in self.hi_prefixes:
                        if px.name in (
                            None,
                            "",
                        ):
                            continue
                        if px.name not in prefixes:
                            continue
                        if junit.prefix_dimension != 1:
                            px **= junit.prefix_dimension
                        if px.value * junit.value * div_lim1 * epsf > lim1:
                            break
                        xscale = divlim
                        lim0_ = px.value * junit.value * div_lim(xscale)
                        if lim0_ * epsf < lim0:
                            continue
                        exponent = False
                        ranges.append(Range(lim0, junit, px, exponent))
                        lim0 = lim0_
                        if lim0 > lim1 * epsf:
                            break
                    if lim0 > lim1 * epsf:
                        break
                lim = lim1
            jpos = len(ranges)
            if lim0 < lim1 * epsf:
                for junit, jscale in scales[iscale + 1 :]:
                    if not junit.output:
                        continue
                    prefixes = junit.out_prefix
                    for px in self.lo_prefixes:
                        if junit.prefix_dimension != 1:
                            px **= junit.prefix_dimension
                        if px.name in (
                            None,
                            "",
                        ):
                            continue
                        if px.name not in prefixes:
                            continue
                        if px.value * junit.value > lim1:
                            continue
                        xscale = divlim
                        lim1_ = px.value * junit.value * div_lim(xscale)
                        if lim1_ < lim0:
                            break
                        lim1 = max(lim0, px.value * junit.value * div_lim1)
                        exponent = False
                        ranges.insert(jpos, Range(lim1, junit, px, exponent))
                        xscale = divlim
                        lim = unit.value * div_lim(xscale) / xscale
                        if lim0 > lim1 * epsf:
                            break
                    if lim0 > lim1 * epsf:
                        break

            # print('\n====')
            # for r in ranges:
            #     print(r)
            # print('x', lim0, lim1)

            # fill in gap with scales
            if lim0 < lim1 * epsf:
                exponent = True
                r = ranges[jpos - 1]
                if r.prefix_value == 1 or r.unit.hi_prefix_power:
                    prefix = r.prefix
                    u = r.unit
                else:
                    prefix = Prefix.neutral
                    u = unit_
                if u.hi_power_limit is not None:
                    hi_digits = int(ceil(log10(u.hi_power_limit * epsm)))
                    lim0_ = _div_lim(u.hi_power_limit, digits - hi_digits)
                    # lim0_ = u.hi_power_limit * div_lim1
                    if lim0_ > lim0:
                        ranges.insert(jpos, Range(lim0, u, prefix, exponent))
                        lim0 = lim0_
                        jpos += 1
                else:
                    ranges.insert(jpos, Range(lim0, u, prefix, exponent))
                    lim0 = lim1
                    jpos += 1

            # print('\n---')
            # for r in ranges:
            #     print(r)
            # print('x', lim0, lim1)

            if lim0 < lim1 * epsf:
                exponent = True
                u = None
                if len(ranges) > jpos:
                    r = ranges[jpos]
                    if r.prefix_value == 1 or r.unit.lo_prefix_power not in (
                        False,
                        "+",
                    ):
                        prefix = r.prefix
                        u = r.unit
                if u is None:
                    prefix = Prefix.neutral
                    u = unit

                if u.lo_power_limit is not None:
                    lim1_ = u.lo_power_limit * div_lim1
                    if lim1_ < lim1:
                        lim1 = max(lim0, lim1_)
                        ranges.insert(jpos, Range(lim1, u, prefix, exponent))
                else:
                    lim1 = lim0
                    ranges.insert(jpos, Range(lim1, u, prefix, exponent))
            if lim0 < lim1 * epsf:
                u = self.units.base
                ranges.insert(jpos, Range(lim0, u, prefix, exponent))
                lim1 = lim0

            prefix = Prefix.neutral
            exponent = False
            ranges.append(Range(lim, unit, prefix, exponent))
            lim = unit.value * div_lim(scale)
            unit_ = unit

        # add high values
        hi_prefix_base = self.hi_prefix_base
        # breakpoint()
        if hi_prefix_base is None:
            for unit in self.outunits[::-1]:
                if len((unit.out_prefix & self.hi_prefixes_set) - set(("",))) > 0:
                    hi_prefix_base = unit
                    break
        if hi_prefix_base is not None:
            unit = hi_prefix_base
            if isinstance(unit, str):
                unit = self.units[unit]
            scale_ = divlim
            for prefix in self.hi_prefixes:
                if unit.prefix_dimension != 1:
                    prefix **= unit.prefix_dimension
                if prefix.name in ("", None):
                    continue
                if prefix.name not in unit.out_prefix:
                    continue
                value = unit.value * prefix.value
                lim = value * div_lim1  # (1, 0) # just to be sure use 0.95
                r = ranges[-1]
                if lim <= r.start:
                    continue
                value_ = r.value
                if scale_ is not None:
                    scale = scale_
                    scale_ = None
                else:
                    scale = value / value_
                if scale > divlim:
                    scale = divlim
                lim = value_ * div_lim(scale)
                exponent = False
                ranges.append(Range(lim, unit, prefix, exponent))

        # just power beyond prefixes
        r = ranges[-1]
        lim = r.start
        unit = r.unit
        prefix = r.prefix
        if self.hi_power_base is not None:
            unit = self.hi_power_base
            prefix = Prefix.neutral
            if isinstance(unit, str):
                try:
                    unit = self.units[unit]
                except:
                    scale, power, prefix, unit = self.unit2sppu(unit)
        exponent = True
        lim *= divlim
        ranges.append(Range(lim, unit, prefix, exponent))

        # add low values
        lo_prefix_base = self.lo_prefix_base
        if lo_prefix_base is None:
            for unit in self.outunits:
                if len((unit.out_prefix & self.lo_prefixes_set) - set(("",))) > 0:
                    lo_prefix_base = unit
                    break
        if lo_prefix_base is not None:
            unit = lo_prefix_base
            if isinstance(unit, str):
                unit = self.units[unit]
            for prefix in self.lo_prefixes:
                if unit.prefix_dimension != 1:
                    prefix **= unit.prefix_dimension
                if prefix.name in ("", None):
                    continue
                if prefix.name not in unit.out_prefix:
                    continue
                value = unit.value * prefix.value
                r = ranges[0]
                lim_ = r.start
                if value > lim_:
                    continue
                value_ = r.value
                scale = value_ / value
                if scale > divlim:
                    scale = divlim
                    # this was added for symmetry with high power, not sure it is correct
                lim_ = value * div_lim(scale)
                r.set_start(lim_)
                lim = (
                    value * div_lim1
                )  # preliminary, maybe overwritten later (see line above)
                exponent = False
                ranges.insert(0, Range(lim, unit, prefix, exponent))

        r = ranges[0]
        unit = r.unit
        prefix = r.prefix
        lim_ = r.value * div_limd
        ranges[0].set_start(lim_)
        if self.lo_power_base is not None:
            unit = self.lo_power_base
            prefix = Prefix.neutral
            if isinstance(unit, str):
                try:
                    unit = self.units[unit]
                except:
                    scale, power, prefix, unit = self.unit2sppu(unit)
        exponent = True
        lim = -inf
        ranges.insert(0, Range(lim, unit, prefix, exponent))

        if self.debug:
            for r in ranges:
                print(r)

        return ranges

    ########################################################################

    def _compute_ranges_v2(self, digits, dec_lim, /):

        divlim = 10 ** (1 + digits)

        def div_lim(scale):
            scale *= dec_lim
            mag = int(ceil(log10(scale)))
            decimals = digits - max(1, mag)  # leading 0 counts as digit
            if mag > digits:
                decimals += 1  # drop decimal point
            decimals += max(0, -int(floor(log10(dec_lim))))
            return _div_lim(scale, decimals)

        div_lim1 = div_lim(1)
        div_lim1d = div_lim(10**digits) / 10**digits
        div_limd1 = div_lim(10 ** (digits + 1))

        # compute divlim scales.
        # set up output scalings as defined
        scales = []
        unit_ = None
        for unit in self.outunits:
            if unit_ is not None:
                if unit_.div_lim is not None:
                    scale = min(unit_.div_lim, divlim)
                else:
                    scale = min(unit.value / unit_.value, divlim)
                    iscale = int(scale)
                    if iscale == scale:
                        scale = iscale
                scales.append((unit_, scale))
            unit_ = unit
        else:
            if unit.div_lim is not None:
                scale = unit.div_lim
            else:
                scale = divlim
            scales.append((unit, scale))

        # set up ranges
        ranges = list()

        # compute ranges for scales
        for xscale in scales:
            unit, scale = xscale
            start_ = unit.value * div_lim1  # smallest value to round to 1
            stop_ = unit.value * div_lim(scale)
            ranges.append(Range(start_, stop_, unit))

        # add phoney boundary ranges
        ranges.insert(0, Range(-inf, 0, Unit(out_prefix=())))
        ranges.append(Range(+inf, +inf, Unit(out_prefix=())))

        # fill in gaps with prefixes
        start = 0
        ranges0 = ranges.copy()
        j = -1
        for i, r in enumerate(ranges0):
            j += 1
            stop = r.start
            if start < stop:
                # from below
                added = 0
                for r0 in ranges0[max(0, i - 1) :: -1]:
                    if not r0.pure:
                        continue
                    unit = r0.unit
                    for prefix in self.hi_prefixes:
                        if prefix.name not in unit.out_prefix:
                            continue
                        if unit.prefix_dimension != 1:
                            prefix **= unit.prefix_dimension
                        start_ = prefix.value * unit.value * div_lim1
                        if start_ >= stop * epsf:
                            break
                        stop_ = prefix.value * unit.value * div_limd1
                        if stop_ <= start * epsm:
                            continue
                        ranges.insert(j, Range(start_, stop_, unit, prefix))
                        j += 1
                        start = stop_
                        added += 1
                    if added > 0:
                        break
                # from above
                added = 0
                for r0 in ranges0[i : i + 1]:  # [i:]
                    if not r0.pure:
                        continue
                    unit = r0.unit
                    for prefix in self.lo_prefixes:
                        if prefix.name not in unit.out_prefix:
                            continue
                        if unit.prefix_dimension != 1:
                            prefix **= unit.prefix_dimension
                        start_ = prefix.value * unit.value * div_lim1
                        if start_ >= stop * epsf:
                            continue
                        stop_ = prefix.value * unit.value * div_limd1
                        if stop_ <= start * epsm:
                            break
                        ranges.insert(j, Range(start_, stop_, unit, prefix))
                        stop = start_
                        added += 1
                    if added > 0:
                        j += added
                        break
            start = r.stop

        # fill hi/lo with special unit/prefixes if requested
        if (start := ranges[0].stop) < (stop := ranges[1].start):
            unit = self.lo_prefix_base
            if isinstance(unit, str):
                unit = self.units[unit]
            if unit is None:
                # TODO - search for lowest unit that allows prefix if base does not
                unit = self.units.base

            for prefix in self.lo_prefixes:
                if prefix.name not in unit.out_prefix:
                    continue
                if unit.prefix_dimension != 1:
                    prefix **= unit.prefix_dimension
                start_ = prefix.value * unit.value * div_lim1
                if start_ >= stop * epsf:
                    continue
                stop_ = prefix.value * unit.value * div_limd1
                if stop_ <= start * epsm:
                    break
                ranges.insert(1, Range(start_, stop_, unit, prefix))
                stop = start_
        if (start := ranges[-2].stop) < (stop := ranges[-1].start):
            unit = self.hi_prefix_base
            if isinstance(unit, str):
                unit = self.units[unit]
            if unit is None:
                # TODO - search for highest unit that allows prefix if base does not
                unit = self.units.base
            for prefix in self.hi_prefixes:
                if prefix.name not in unit.out_prefix:
                    continue
                if unit.prefix_dimension != 1:
                    prefix **= unit.prefix_dimension
                start_ = prefix.value * unit.value * div_lim1
                if start_ >= stop * epsf:
                    break
                stop_ = prefix.value * unit.value * div_limd1
                if stop_ <= start * epsm:
                    continue
                ranges.insert(1, Range(start_, stop_, unit, prefix))
                start = stop_

        # fill in with (prefix) powers of directly(!) neighboring ranges
        start = 0
        ranges0 = ranges.copy()
        j = -1
        for i, r in enumerate(ranges0):
            j += 1
            stop = r.start
            # from below
            while start < stop:
                r0 = ranges0[i - 1]
                if r0.stop <= 0:
                    break
                unit = r0.unit
                prefix = r0.prefix
                stop_ = unit.hi_power_limit
                if stop_ is None:
                    stop_ = +inf
                else:
                    # this needs to be more sophisticated;
                    # add extra flag to range
                    stop_ *= div_lim1d
                if stop_ < start:
                    break
                if r0.has_prefix and not unit.hi_prefix_power:
                    break
                ranges.insert(j, Range(start, stop_, unit, prefix, exponent=True))
                j += 1
                start = stop_
                break
            # from above
            while start < stop:
                r0 = ranges0[i]
                if r0.start >= inf:
                    break
                unit = r0.unit
                if r0.has_prefix and not unit.lo_prefix_power:
                    break
                prefix = r0.prefix
                start_ = unit.lo_power_limit
                if start_ is None:
                    start_ = 0
                else:
                    # this needs to be more sophisticated;
                    # add extra flag to range
                    start_ *= div_lim1
                if start_ > stop:
                    break
                ranges.insert(j, Range(start_, inf, unit, prefix, exponent=True))
                j += 1
                stop = start_
                break

            start = r.stop

        # add hi/lo ranges (if needed)
        # TODO add xx_power_limit treatment
        if (start := ranges[0].stop) < ranges[1].start:
            unit = self.lo_power_base
            if unit is None:
                unit = self.units.base
            prefix = Prefix.neutral
            if isinstance(unit, str):
                try:
                    unit = self.units[unit]
                except:
                    scale, power, prefix, unit = self.unit2sppu(unit)
            ranges.insert(1, Range(start, inf, unit, prefix, True))
        if ranges[-1].start > (start := ranges[-2].stop):
            unit = self.hi_power_base
            if unit is None:
                unit = self.units.base
            prefix = Prefix.neutral
            if isinstance(unit, str):
                try:
                    unit = self.units[unit]
                except:
                    scale, power, prefix, unit = self.unit2sppu(unit)
            ranges.insert(-1, Range(start, inf, unit, prefix, True))

        # fill remaining gaps with default unit
        start = 0
        ranges0 = ranges.copy()
        j = -1
        unit = self.units.base
        prefix = Prefix.neutral
        for i, r in enumerate(ranges0):
            j += 1
            stop = r.start
            if start < stop:
                ranges.insert(j, Range(start, inf, unit, prefix, True))
                j += 1
            start = r.stop

        # clean overlaps
        # TODO!!!

        if self.debug:
            for r in ranges:
                print(r)
            print(len(ranges))

        return ranges[1:-1]

    def unit2sppu(self, unit):
        if unit in (
            ...,
            "",
        ):
            unit = self.units.base
        if isinstance(unit, Unit):
            unit = unit.name
        value, power, unit = self.split_unit3(unit)
        assert float(value) == 1, f'Require pure power as unit but received "{value}".'
        power = Power(value=power)
        scale, unit, prefix = self._unit2scale.get(unit)
        scale *= power.value
        return scale, power, prefix, unit

    def __call__(
        self,
        value,
        cut=True,
        base=None,
        strip=True,
        dec_lim=None,
        unicode=False,
        latex=False,
        comma=False,
        numeric_int=False,
        rounding=False,
        digits=None,
        sign=False,
        align=False,
        return_string=True,
        return_range=False,
        return_mantissa=False,
        return_mantissa_numeric=False,
        return_unit=False,
        return_value=False,
        unit=None,
        unit_upgrade=False,  # allow switch to larger scale
    ):
        """
        more parameters should set defaults in __init__
        """

        if digits is None:
            digits = self.digits
        if dec_lim is None:
            dec_lim = self.dec_lim

        ranges_key = (digits, dec_lim)
        ranges = self.ranges.get(ranges_key, None)
        if ranges is None:
            ranges = self._compute_ranges(digits, dec_lim)
            self.ranges[ranges_key] = ranges

        val = value
        if isinstance(val, str):
            val = self.human2val(val)

        if base is None:
            base = self.base

        if self.convert is not None:
            val = self.convert(val)
        elif isinstance(base, str):
            val *= self.unit2scale(base)
        elif base is not None:
            val *= base
        else:
            val = val

        aval = abs(val)
        if self.integer:
            oval = round(aval)
            assert abs((oval + 1 / (2 << 53)) / (aval + 1 / (2 << 53)) - 1) < 1 / (
                2 << 46
            )
            xval = oval
        else:
            xval = aval

        if unit is not None:
            scale, power, prefix, unit = self.unit2sppu(unit)
            ru = Range(-1, unit, prefix, False, power)
            if unit_upgrade is True:
                divlim = _div_lim(10 ** (digits + 1))
                r = ranges[0]
                for r_ in ranges[1:]:
                    if xval <= r.value * divlim or xval <= r_.start:
                        break
                    r = r_
                if r.start <= scale:
                    r = ru
            elif unit_upgrade == "range":
                # search duplicate copy from below, but may be modified
                r = ranges[0]
                for r_ in ranges[1:]:
                    if xval < r_.start:
                        break
                    r = r_
                if r.start <= scale:
                    r = ru
            else:
                r = ru

        if unit is None:
            # TODO - add faster search
            r = ranges[0]
            for r_ in ranges[1:]:
                if xval < r_.start:
                    break
                r = r_

        xval /= r.value

        if xval == 0 and unit is None:
            scale, power, prefix, unit = self.unit2sppu(self.zero_unit)
            r = Range(-1, unit, prefix, False, power)

        if xval > 0 and rounding:
            xval = round(xval, digits - 1 - floor(log10(xval)))

        if xval == 0:
            value = Value.zero
        elif r.exponent:
            sv = f"{xval:{digits+5}.{digits-1}e}"
            m, e = sv.split("e")
            e = int(e)
            r += Power(e)
            value = Value(m)
        else:
            value = Value(xval)

        m = value.str(
            digits=digits,
            comma=comma,
            cut=cut,
        )

        if val < 0:
            m = "-" + m
        elif sign is True:
            m = "+" + m

        if align:
            m = f"{m:>{digits+2}s}"

        s = m + r.str(unicode=unicode, latex=latex)

        if latex and s.startswith(r"1\times"):
            s = s[len(r"1\times") :]

        retval = list()
        if return_string:
            retval.append(s)
        if return_range:
            retval.append(r)
        if return_mantissa:
            retval.append(m)
        if return_mantissa_numeric:
            numeric = xval
            if val < 0:
                numeric = -numeric
            if numeric_int:
                if m.find(".") == -1:
                    numeric = int(numeric)
            retval.append(numeric)
        if return_unit:
            retval.append(r.str(unicode=unicode, latex=latex, unit=True))
        if return_value:
            retval.append(r.value)

        if len(retval) == 1:
            return retval[0]
        return tuple(retval)

    def all_units(self):
        if not hasattr(self, "_all_units"):
            all_units = []
            for units in self._unit2scale.keys():
                for special in "|*.()[]{}?^$":
                    units = units.replace(special, "\\" + special)
                all_units.append(units)
            self._all_units = "|".join(all_units)
        return self._all_units

    def split_pattern(self, s):
        if not hasattr(self, "_split_pattern"):
            pattern = rf"^\s*([0-9.]+)?([EeDd][+-]?[0-9]+)?\s*({self.all_units()})?\s*$"
            self._split_pattern = re.compile(pattern)

        return self._split_pattern.findall(s)[0]

    def split_unit3(
        self,
        s,
        num_val=False,
        num_unit=False,
        num_power=False,
    ):
        if not hasattr(self, "_split_pattern"):
            pattern = rf"^\s*([0-9.]+)?([EeDd][+-]?[0-9]+)?\s*({self.all_units()})?\s*$"
            self._split_pattern = re.compile(pattern)

        value, power, unit = self.split_pattern(s)
        if value == "":
            value = "1"
        if value.count(".") > 0:
            value = value.rstrip("0").rstrip(".")
        if power != "":
            power = power.replace("d", "e").replace("D", "E")
        if num_val:
            try:
                value = int(value)
            except:
                value = float(value)
            ivalue = int(value)
            if value == ivalue:
                value = ivalue
        if num_unit:
            unit = self.unit2scale(unit)
        if num_power:
            power = Power(value=power).power
        return value, power, unit

    def split_unit2(self, s, num_val=False, num_unit=False):
        if not hasattr(self, "_split_pattern"):
            pattern = rf"^\s*([0-9.]+)?([EeDd][+-]?[0-9]+)?\s*({self.all_units()})?\s*$"
            self._split_pattern = re.compile(pattern)

        value, power, unit = self.split_pattern(s)
        if value == "":
            value = "1"
        if value.count(".") > 0:
            value = value.rstrip("0").rstrip(".")
        if power != "":
            value += power.replace("d", "e").replace("D", "E")
        if num_val:
            try:
                value = int(value)
            except:
                value = float(value)
            ivalue = int(value)
            if value == ivalue:
                value = ivalue
        if num_unit:
            unit = self.unit2scale(unit)
        return value, unit

    def split_unit(self, s, num_val=False, num_unit=False):
        value, power, unit = self._split_pattern.findall(s)[0]
        if value == "":
            value = "1"
        if num_val:
            try:
                value = int(value)
            except:
                value = float(value)
            ivalue = int(value)
            if value == ivalue:
                value = ivalue
        else:
            if value.count(".") > 0:
                value = value.rstrip("0").rstrip(".")
        if num_unit:
            unit = self.unit2scale(unit)
        if power != "":
            power = "1" + power.replace("d", "e").replace("D", "E")
            if num_unit:
                unit *= float(power)
            else:
                unit = f"{power} {unit}"
        return value, unit

    def unit2scale(self, unit, base=None):
        try:
            power = float(unit)
            unit = ""
        except ValueError:
            power = None

        if base is None:
            base = self.base

        if power is None:
            if not hasattr(self, "_unit_pattern"):
                pattern = (
                    rf"^\s*([0-9.]+(?:[EeDd][+-]?[0-9]+)?)?\s*({self.all_units()})?\s*$"
                )
                self._unit_pattern = re.compile(pattern)
            try:
                power, unit = self._unit_pattern.findall(unit)[0]
            except IndexError:
                raise AttributeError(f'Cannot understand unit "{unit}".')
        if unit != "":
            scale = self._unit2scale.get(unit, (None,))[0]
            if scale is None:
                raise AttributeError(f'Cannot understand unit "{unit}".')
        else:
            scale = 1
        if isinstance(base, str):
            scale /= self.unit2scale(base)
        elif base is not None:
            scale /= base
        if scale != 0:
            if power != "":
                scale *= float(power)
            iscale = int(scale)
            if scale == iscale:
                scale = iscale
        return scale

    def max_unit(self, units):
        return sorted(units, key=self.unit2scale)[-1]

    def human2val(self, human, base=None):
        if not isinstance(human, str):
            if isinstance(human, int):
                value = human
            else:
                value = float(human)
            unit = self.zero_unit.name
        else:
            value, unit = self.split_unit2(human)
            if unit == "":
                unit = self.units.base.name
        scale = self.unit2scale(unit, base)
        value = float(value)
        value *= scale
        if value < (2 << 53):
            ivalue = int(value)
            if value == ivalue:
                value = ivalue
        return value

    def __getattr__(self, attr):
        try:
            return self._unit2scale[attr][0]
        except KeyError:
            raise AttributeError
