import numpy as np


def version2human(version, strip=True, digits=2):
    """
    Convert version number to string.

    20000 becomes " 2.00.00"

    PARAMETERS
    strip [True]: remove trailing space
    """
    vr = version
    scale = 10**digits
    scale2 = scale**2
    if isinstance(vr, np.float64):
        vr = int(round(vr * scale2))
    v1 = vr // scale2
    vr = vr - scale2 * v1
    v2 = vr // scale
    v3 = vr - v2 * scale
    s = "{1:{0}d}.{2:0{0}d}.{3:0{0}d}".format(digits, v1, v2, v3)
    if strip:
        s = s.strip()
    return s
