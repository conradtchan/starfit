#! /bin/env python3

import datetime
import sys

try:
    import physconst
except ImportError:
    SEC = 31556926
else:
    SEC = physconst.SEC

_Units = ("", "k", "M", "G", "T", "P", "E", "Z", "Y")
_units = ("", "m", "u", "n", "p", "f", "a")


def _div_lim(x, digits=0):
    return x * (1 - 2.0e-15) - 0.5 * 10 ** (-digits)


def time2human(time, digits=2, cut=True, extended=False, strip=True):
    """Convert time in seconds in readable format.

    Specify number of total *digits*
    and whether to *cut* trailing zeros.

    If "extended" is set to True, also the numeric value,
    the unit string and the scale factor are returned.
    """

    if isinstance(time, datetime.timedelta):
        time = time.total_seconds()

    atime = abs(time)
    xtime = atime
    su = "s"

    div_lim1000 = _div_lim(1000)
    div_lim100 = _div_lim(100)
    div_lim60 = _div_lim(60)
    div_lim1 = _div_lim(1, 3)

    if atime > 1:
        # big numbers
        if atime > div_lim100:
            xtime /= 60
            su = "min"
            if atime > 60 * div_lim60:
                xtime /= 60
                su = "h"
                if atime > 24 * 60 * div_lim100:
                    xtime /= 24
                    su = "d"
                    if atime > SEC:
                        xtime = atime / SEC
                        su = "yr"
                        i = 0
                        while xtime > div_lim1000:
                            xtime /= 1000
                            i += 1
                        if i >= len(_Units):
                            return "***"
                        su = _Units[i] + su
    elif atime > 0:
        # small numbers ...
        i = 0
        while xtime < div_lim1:
            xtime *= 1000
            i += 1
        if i >= len(_units):
            return "***"
        su = _units[i] + su

    sv = f"{xtime:20.15f}".strip()
    i = sv.find(".")
    l = max(digits + 1, i + 1)
    format = "{:" + f"{l:d}".strip() + "." + f"{l - i - 1:d}".strip() + "f}"
    sv = format.format(xtime).strip()

    if cut:
        if sv.find(".") > 0:
            sv = sv.rstrip("0")
        sv = sv.rstrip(".")
    if time < 0:
        sv = "-" + sv
    s = sv + " " + su

    if strip:
        s = s.strip()
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


if __name__ == "__main__":
    argv = sys.argv
    if len(argv) == 2:
        try:
            print(time2human(float(argv[1])))
        except:
            print("***")
