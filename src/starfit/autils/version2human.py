import numpy as np


def version2human(version, strip=True):
    """
    Convert version number to string.

    20000 becomes " 2.00.00"

    PARAMETERS
    strip [True]: remove trailing space
    """
    vr = version
    if isinstance(vr, np.float64):
        vr = int(round(vr * 10000))
    v1 = vr // 10000
    vr = vr - 10000 * v1
    v2 = vr // 100
    v3 = vr - v2 * 100
    s = f"{v1:2d}.{v2:02d}.{v3:02d}"
    if strip:
        s = s.strip()
    return s
