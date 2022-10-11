#! /bin/env python3


# TODO - check consistent use of digits for div_lim
# TODO - check use of div_lim for rounding
# TODO - check use of div_lim, dec_lim, rounding for unit_upgrade

import datetime
import math
import sys

from .util import _div_lim, _prefixes, _Prefixes

try:
    from ..physconst import YR
except ImportError:
    YR = 31556926

_unit_base = dict(
    s=1,
    min=60,
    h=3600,
    d=86400,
    yr=YR,
)


def time2human(
    time,
    digits=2,
    cut=True,
    extended=False,
    strip=True,
    unit=None,
    unit_upgrade=False,
    rounding=False,
    comma=False,
    numeric_int=False,
    dec_lim=None,
    fixed=False,
    base=None,
):
    # dec_lim = 10 results in integer output
    """Convert time in seconds in readable format.

    Specify number of total *digits*
    and whether to *cut* trailing zeros.

    If "extended" is set to True, also the numeric value, the unit string, and the
    scale factor are returned, and the numeric value can be returned
    as integer (numeric_int).

    A minimum decimal limit (dec_lim) can be set.  With default
    setting, if this is set to 10, no decimals are produced.

    The output *unit* can be enforced as well as *rounding* to the
    specified number of digits and adding of *comma*s.

    Experimentally, *base* is the unit of the input (default: seconds)
    """

    if isinstance(time, datetime.timedelta):
        time = time.total_seconds()
    elif isinstance(base, str):
        time *= unit2scale(base)
    elif base is not None:
        time *= base

    atime = abs(time)
    xtime = atime
    su = "s"

    div_lim1000 = _div_lim(1000)
    div_lim100 = _div_lim(100)
    div_lim60 = _div_lim(60)
    div_lim1 = _div_lim(1, 3)

    if dec_lim is None:
        div_lims = [div_lim1, div_lim100, div_lim60, div_lim100, div_lim1000]
        div_lim = 1
    else:
        div_lims = [
            _div_lim(dec_lim, 3),
            _div_lim(dec_lim * 60),
            _div_lim(dec_lim * 60),
            _div_lim(dec_lim * 100),
            _div_lim(dec_lim * 1000),
        ]
        div_lim = dec_lim

    decimals = 0
    if unit is not None:
        su = unit
        xtime /= unit2scale(su)
        if unit_upgrade:
            while True:
                if (su[0] in _prefixes[1:]) and (su != "min") and (xtime > div_lim1000):
                    i = _prefixes.index(su[0])
                    su = _prefixes[i - 1] + su[1:]
                    xtime /= 1000
                    decimals += 3
                else:
                    break
            if su in _unit_base and (xtime > div_lim1000):
                su = _Prefixes[1] + su
                xtime /= 1000
                decimals += 3
            while True:
                if (su[0] in _Prefixes[1:]) and (xtime > div_lim1000):
                    i = _Prefixes.index(su[0])
                    su = _Prefixes[i + 1] + su[1:]
                    xtime /= 1000
                    decimals += 3
                else:
                    break
    elif atime >= div_lim:
        # big numbers
        if atime > div_lims[1]:
            xtime /= 60
            su = "min"
            if atime > 60 * div_lims[2]:
                xtime /= 60
                su = "h"
                if atime > 24 * 60 * div_lims[3]:
                    xtime /= 24
                    su = "d"
                    if atime > YR * div_lims[0]:
                        xtime = atime / YR
                        su = "yr"
                        i = 0
                        while xtime > div_lims[4]:
                            xtime /= 1000
                            i += 1
                        if i >= len(_Prefixes):
                            return "***"
                        su = _Prefixes[i] + su
    elif atime > 0:
        # small numbers ...
        i = 0
        while xtime < div_lims[0]:
            xtime *= 1000
            i += 1
        if i >= len(_prefixes):
            return "***"
        su = _prefixes[i] + su

    sv = f"{xtime:20.15f}".strip()
    i = sv.find(".")
    l = max(digits + 1, i + 1) + decimals
    format = "{:" + f"{l:d}".strip() + "." + f"{l - i - 1:d}".strip() + "f}"

    if xtime > 0 and rounding:
        xtime = round(xtime, digits - 1 - math.floor(math.log10(xtime)))

    sv = format.format(xtime).strip()

    if cut:
        if sv.find(".") > 0:
            sv = sv.rstrip("0")
        sv = sv.rstrip(".")

    if comma:
        l = len(sv)
        j = sv.find(".")
        if j == -1:
            j = l
        for i in range(j - 3, 0, -3):
            sv = sv[:i] + "," + sv[i:]
    if time < 0:
        sv = "-" + sv
    if strip:
        sv = sv.strip()
    if fixed:
        if len(sv) < digits + 1:
            sv = f"{sv:>{digits + 1}}"
    s = sv + " " + su
    if extended:
        if xtime == 0:
            scale = 1
        else:
            scale = atime / xtime
        numeric = time / scale
        #        if sv.find('.') == -1:
        #            numeric = int(numeric)
        return s, numeric, su, scale
    return s


def split_unit(s, num_val=False, num_unit=False):
    if s.count(" ") == 1:
        v, u = s.split()
    else:
        j = -1
        for i, c in enumerate(s):
            if c in "1234567890.":
                j = i + 1
        if j == -1:
            v = 1
            u = s.strip()
        elif j == len(s):
            v = s.strip()
            u = "s"
        else:
            v = s[:j].strip()
            u = s[j:].strip()
    if num_val:
        try:
            v = int(v)
        except:
            v = float(v)
        iv = int(v)
        if iv == v:
            v = iv
    if num_unit:
        u = unit2scale(u)
    return v, u


_unit2scale = dict()
for b, v in _unit_base.items():
    for i, u in enumerate(_prefixes):
        _unit2scale[u + b] = v * 10 ** (-3 * i)
    for i, u in enumerate(_Prefixes[1:]):
        _unit2scale[u + b] = v * 10 ** (3 * (i + 1))


def unit2scale(unit, base=None):
    scale = _unit2scale.get(unit, 0)
    if isinstance(base, str):
        scale /= unit2scale(base)
    elif base is not None:
        scale /= base
    if scale != 0:
        return scale


def max_unit(units):
    return sorted(units, key=unit2scale)[-1]


def human2time(s, base=None):
    s = s.replace(",", "")  # ???
    v, u = split_unit(s)
    try:
        return int(v) * unit2scale(u, base)
    except:
        return float(v) * unit2scale(u, base)


if __name__ == "__main__":
    argv = sys.argv
    if len(argv) == 2:
        try:
            print(time2human(float(argv[1])))
        except:
            print("***")
    elif len(argv) == 3:
        try:
            base = float(argv[2])
        except:
            base = argv[2]
        try:
            print(time2human(float(argv[1]), base=base))
        except:
            print("***")
