"""
some helper tools and common definitions

Unofficial:
   zepto (1e-21)
   yocto (1e-24)
"""

_Prefixes = ("", "k", "M", "G", "T", "P", "E", "Z", "Y")
_prefixes = ("", "m", "u", "n", "p", "f", "a", "z", "y")
_prefixes_unicode = ("", "m", "\u03BC", "n", "p", "f", "a", "z", "y")
_prefixes_latex = ("", "m", r"\mu", "n", "p", "f", "a", "z", "y")

_BinPrefixes = ("", "ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi")
_binPrefixes = ("", "mi", "ui", "ni", "pi", "fi", "ai", "zi", "yi")
_binPrefixes_unicode = ("", "mi", "\u03BCi", "ni", "pi", "fi", "ai", "zi", "yi")
_binPrefixes_latex = ("", "mi", r"{\mu}i", "ni", "pi", "fi", "ai", "zi", "yi")

_BinPrefixes_power = 1024


# this may be good for float64 only
def _div_lim(x, digits=0):
    return x * (1 - 2 ** (-53)) - 0.5 * 10 ** (-digits)
