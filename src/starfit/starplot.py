"""Various plotting routines"""

from distutils.spawn import find_executable

import matplotlib as mpl

found_latex = find_executable("latex")
found_dvipng = find_executable("dvipng")

if found_latex and found_dvipng:
    mpl.rc("text", usetex=True)
else:
    mpl.rc("font", family="serif")
    mpl.rc("font", serif=["DejaVu Serif", "Computer Modern Roman"])
    mpl.rc("mathtext", fontset="cm")

if not found_latex:
    print(
        "\nWarning: LaTeX installation not found."
        "Reverting to use of 'mathtext' for plots.\n"
    )
if not found_dvipng:
    print(
        "\nWarning: dvipng installation not found."
        "Reverting to use of 'mathtext' for plots.\n"
    )

import numpy as np

from .autils.keputils import mass_string

_plot_unit_translate = {
    "solar masses": r"$\mathsf{M}_\odot$",
    "M_sun": r"$\mathsf{M}_\odot$",
    "Msun": r"$\mathsf{M}_\odot$",
    "Z_sun": r"$\mathsf{Z}_\odot$",
    "Zsun": r"$\mathsf{Z}_\odot$",
    "solar": r"$\mathsf{Z}_\odot$",
    "He core fraction": r"$\mathsf{M}_{\mathsf{He}}$",
}


def _plot_unit_format_mixing(value, *_):
    if value == 0:
        return "(no mixing)"
    return rf"$\log(f_{{\mathsf{{mix}}}})={np.log10(value):4.1f}$"


_plot_unit_formatters = {
    "model": lambda value, *_: f"Model {value}",
    "He core fraction": _plot_unit_format_mixing,
    "solar masses": lambda value, *_: rf"{mass_string(value)} $\mathsf{{M}}_{{\odot}}$",
}

_plot_title_formatters = {
    "lower mass cut": lambda value, *_: rf"$\ge{mass_string(value)}\,\mathsf{{M}}_{{\odot}}$",
    "upper mass cut": lambda value, *_: rf"$\le{mass_string(value)}\,\mathsf{{M}}_{{\odot}}$",
    "gamma": lambda value, form, *_: rf"$\Gamma={value:{form}}$",
    "eexp": lambda value, form, *_: rf"$E^{{\mathsf{{exp}}}}={value:{form}}$",
    "sigma": lambda value, form, unit, *_: rf"$\pm{value:{form}}$ {unit}",
}


def leg_starname(ax, name):
    legt = ax.legend(
        [],
        [],
        title=name,
        title_fontproperties=dict(
            size="x-large",
        ),
        bbox_to_anchor=[0.99, 0.99],
        loc="upper right",
        frameon=False,
    )
    legt.set_draggable(True)
    ax.add_artist(legt)


def leg_copyright(ax):
    legc = ax.legend(
        [],
        [],
        title=r"$\copyright$ www.StarFit.org",
        title_fontproperties=dict(
            size="x-small",
        ),
        bbox_to_anchor=[-0.005, -0.02],
        loc="lower left",
        frameon=False,
    )
    legc.set_draggable(True)
    ax.add_artist(legc)


def leg_info(ax, lines):
    if isinstance(lines, str):
        lines = lines.splitlines()
    lines = "\n".join(lines)
    legc = ax.legend(
        [],
        [],
        title=lines,
        title_fontproperties=dict(
            size="medium",
        ),
        bbox_to_anchor=[0.99, 0.05],
        loc="lower right",
        frameon=False,
    )
    legc.set_draggable(True)
    ax.add_artist(legc)


class IntFormatter(mpl.ticker.FuncFormatter):
    def __init__(self, *args, **kwargs):
        def formatter(*args, **kwargs):
            """
            function to format string for use with ticker
            """
            v = args[0]
            if int(v) == v:
                v = int(v)
                return str(v)
            return f"{v:.2g}"

        super().__init__(formatter)
