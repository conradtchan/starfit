#! /bin/env python3

"""
Module for huma-readable byte strings.
"""

import sys

_units = ("", "k", "M", "G", "T", "P", "E", "Z", "Y")


def byte2human(size, strip=True, SI=False, short=False, length=3):
    """
    Return string in human-readible format.
    """
    assert 2 < length < 6, "Length out of Range"
    su = "B"
    if SI:
        prefix = ""
        div = 1000
    else:
        prefix = "i"
        div = 1024
    div_lim = 1000 * (1 - 0.5 * 10 ** (-length) - 2.0e-15)
    if short:
        prefix = ""
        su = ""

    asize = abs(size)
    osize = round(asize)
    assert abs((osize + 1.0e-16) / (asize + 1.0e-16) - 1) < 1.0e-14
    xsize = osize

    i = 0
    while xsize > div_lim:
        xsize /= div
        i += 1
    if length == 3:
        rnd_lim = 10 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        if (xsize >= rnd_lim) or (i == 0):
            sv = f"{int(round(xsize)):3d}"
        else:
            sv = f"{xsize:3.1f}"
    elif length == 4:
        rnd_lim1 = 100 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        rnd_lim2 = 10 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        if (xsize >= rnd_lim1) or (i == 0):
            sv = f"{int(round(xsize)):4d}"
        elif xsize >= rnd_lim2:
            sv = f"{xsize:4.1f}"
        else:
            sv = f"{xsize:4.2f}"
    elif length == 5:
        rnd_lim1 = 1000 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        rnd_lim2 = 100 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        rnd_lim3 = 10 * (1 - 0.5 * 10 ** (-length + 1) - 2.0e-15)
        if i > 0 and 999 < round(xsize * div) < div:
            xsize *= div
            i -= 1
            sv = f"{int(round(xsize)):4d}"
            sv = sv[0] + "," + sv[1:4]
        elif (xsize >= rnd_lim1) or (i == 0):
            sv = f"{int(round(xsize)):5d}"
        elif xsize >= rnd_lim2:
            sv = f"{xsize:5.1f}"
        elif xsize >= rnd_lim3:
            sv = f"{xsize:5.2f}"
        else:
            sv = f"{xsize:5.3f}"
    else:
        raise Exception("Length out of Range")

    if i >= len(_units):
        sv = "*" * length
    unit = _units[i]
    if i >= 1:
        unit += prefix
    unit = unit + su

    if round(size) < 0:
        sv = "-" + sv.strip()
        sv = " " * (length - len(sv)) + sv
    s = sv + " " + unit
    if strip:
        s = s.strip()

    return s


if __name__ == "__main__":
    argv = sys.argv
    if len(argv) == 2:
        try:
            print(byte2human(round(float(argv[1]))))
        except:
            print("***")
