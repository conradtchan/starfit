import math

# SI DEFINED CONSTANTS
dnucs133hfs = 9192631770  # ?: (genauer ??(133Cs)hfs) Casiumfrequenz, exakt 9192631770 Hz, die Frequenz der Strahlung beim ?bergang zwischen den beiden Hyperfeinstrukturniveaus des Grundzustandes von Atomen des Nuklids 133Cs.
CLIGHT = 29979245800  # : Lichtgeschwindigkeit im Vakuum, exakt 299792458 Meter/Sekunde
PLANCK = 6.62607015e-27  # h: Planck'sches Wirkungsquantum, exakt 6,62607015 * 10-34 Joule*Sekunde
ECHARGE_SI = 1.602176634e-19
ECHARGE = (
    ECHARGE_SI * CLIGHT / 10
)  # e: Elementarladung, exakt 1,602176634 * 10-19 Coulomb,
KB = 1.380649e-16  # k : Boltzmann-Konstante, exakt 1,380649 * 10-23 Joule/Kelvin,
NA = 6.02214076e23  # NA: Avogadro-Konstante, exakt 6,02214076 * 10+23 1/mol
KCD = 683  # Kcd: Lichtausbeute, exakt 683 Lumen/Watt bei monochromatischer Strahlung mit 540 * 1012 Hertz


PI = 4 * math.atan(1.0)
# XMSUN   = 1.9891e33
MSUN = 1.9885e33  # 20220915
LSUN = 3.84e33
# RSUN    = 6.98e10
RSUN = 6.96342e10  # 20220915
# GRAV    = 6.67259e-8
GRAV = 6.67430e-8  # 20220915
PC = 3.08567756706e18
AU = 1.4959787e13
ME = 9.1093897e-28

YR = 31558149.7635456  # siderical year J2000
YR = 31558432.5504  # Anomalistic Year (perihelion) J2000
YR = 31556925.445  # tropical year in the year J2000
YR = 31556952.0  # Gregorian Year
YR = 31556926.0  # tropical Year (old)

# compatibility
XMSUN = MSUN
XLSUN = LSUN
SEC = YR

# SB      = 5.67051e-5
SB = 2 * PI**5 * (KB / PLANCK) ** 3 * KB / (15 * CLIGHT**2)
# KB      = 1.380658e-16
# PLANCK  = 6.6260755e-27
# HBAR    = 1.05457266e-27
HBAR = PLANCK / (2 * PI)
# CLIGHT  = 29979245888e0
# ECHARGE = 4.8032067993e-10
# EV      = 1.60217733e-12
EV = ECHARGE_SI * 1.0e7
KEV = EV * 1.0e3
MEV = EV * 1.0e6
GEV = EV * 1.0e9
# RK      = 83145112.1195e0
RK = KB * NA
# ARAD    = 7.56591414985e-15
ARAD = 4 * SB / CLIGHT
# NA      = 6.0221367e23
# AMU     = 1.6605402e-24
AMU = 1 / NA

DAY = 86400.0
BARN = 1.0e-24
SIGT = 8 * PI / 3 * (ECHARGE**2 / (ME * CLIGHT**2)) ** 2

# masses in MeV
M_N = 939.56542052
M_P = 938.27208816
M_E = 0.5109989500


class Kepler:
    """
    Physics constants as used with KEPLER
    """

    solmass = 1.9892e33
    gee = 6.670e-8
    n0 = 6.02254e23
    sigt = 6.65205e-25
    k = 1.38054e-16
    a = 7.5648e-15
    me = 9.10908e-28
    h = 6.62559e-27
    c = 2.997925e10
    pie = 3.1415926536e0
    solmassi = 1 / solmass
    solrad = 6.9599e10
    solradi = 1 / solrad
    penmex = 0.782333098e0  # (m_n-m_p-m_e)*c**2 [MeV]
    year = 3.1558e7
    pie43 = pie * 4 / 3

    rk = k * n0
    sb = a * c / 4
