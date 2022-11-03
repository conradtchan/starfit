"""Various plotting routines"""

from distutils.spawn import find_executable
from itertools import cycle

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

import colorsys

import matplotlib.pyplot as plt
import numpy as np

from .autils.isotope import ion as I
from .autils.keputils import mass_string
from .autils.stardb import StarDB

_unit_translate = {
    "solar masses": r"$\mathsf{M}_\odot$",
    "M_sun": r"$\mathsf{M}_\odot$",
    "Msun": r"$\mathsf{M}_\odot$",
    "Z_sun": r"$\mathsf{Z}_\odot$",
    "Zsun": r"$\mathsf{Z}_\odot$",
    "solar": r"$\mathsf{Z}_\odot$",
    "He core fraction": r"$\mathsf{M}_{\mathsf{He}}$",
}


def _unit_format_mixing(value, *_):
    if value == 0:
        return "(no mixing)"
    return rf"$\log(f_{{\mathsf{{mix}}}})={np.log10(value):4.1f}$"


_unit_formatters = {
    "model": lambda value, *_: f"Model {value}",
    "He core fraction": _unit_format_mixing,
    "solar masses": lambda value, *_: rf"{mass_string(value)} $\mathsf{{M}}_{{\odot}}$",
}

_title_formatters = {
    "lower mass cut": lambda value, *_: rf"$\ge{mass_string(value)}\,\mathsf{{M}}_{{\odot}}$",
    "upper mass cut": lambda value, *_: rf"$\le{mass_string(value)}\,\mathsf{{M}}_{{\odot}}$",
    "gamma": lambda value, form, *_: rf"$\Gamma={value:{form}}$",
    "eexp": lambda value, form, *_: rf"$E^{{\mathsf{{exp}}}}={value:{form}}$",
    "sigma": lambda value, form, unit, *_: rf"$\pm{value:{form}}$ {unit}",
}


def abuplot(
    indices=None,
    offsets=None,
    star=None,
    database=None,
    database_idx=None,
    database_off=None,
    full_abudata=None,
    eval_data=None,
    list_db=None,
    list_comb=None,
    sun_full=None,
    sun_star=None,
    combine=None,
    database_label=None,
    exclude_index=None,
    uplim_index_star=None,
    lolim_index_all=None,
    solution_fitness=None,
    save=None,
    figsize=(10, 6),
    xlim=None,
    ylim=None,
    show_copyright=True,
    savename=None,
    fontsize=None,
    annosize=None,
    data_size=3,
    dist=None,
):

    zlist_db = np.array([ion.Z for ion in list_db])
    zlist_comb = np.array([ion.Z for ion in list_comb])

    fig = plt.figure(
        figsize=figsize,
        dpi=102,
        facecolor="white",
        edgecolor="white",
    )

    if annosize is None:
        if fontsize is not None:
            annosize = fontsize
        else:
            annosize = "small"

    if fontsize is None:
        fontsize = 12

    ax = fig.add_subplot(111)
    # ax.set_xlabel('Z')
    ax.set_xlabel("Element Charge Number", fontsize=fontsize)
    # ax.set_ylabel('[X]')
    ax.set_ylabel("Logarithm of Abundance Relative to Sun", fontsize=fontsize)

    # First, sum up the values for all the stars in the solution
    if type(indices) is list:
        summed = np.sum((full_abudata[:, indices] + 1.0e-99) * offsets, axis=1)
    else:
        summed = (full_abudata[:, indices] + 1.0e-99) * offsets

    # Transform the list of matched elements to an index for the db
    index_t = np.in1d(zlist_db, zlist_comb, assume_unique=True)

    logsun_full = np.log10(sun_full)
    logsun_star = np.log10(sun_star)

    # The star data points
    y_star = eval_data.abundance - logsun_star
    y_star_det = eval_data.detection - logsun_star

    x_star = np.array([ion.Z for ion in eval_data.element])

    y_star_uco = np.abs(eval_data.error)
    y_star_cov = np.sqrt(np.sum(eval_data.covariance**2, axis=1))
    y_star_err = np.sqrt(y_star_cov**2 + y_star_uco**2)

    y_star_err = np.tile(y_star_err, (2, 1))
    y_star_uco = np.tile(y_star_uco, (2, 1))

    # set asymmetric error for upper limits
    up_lims = uplim_index_star
    y_star_err[1, up_lims] = 0

    if database_label is None:
        database_label = [f"{i:d}" for i in range(len(database))]

    legt = ax.legend(
        [],
        [],
        title=star.name,
        title_fontproperties=dict(
            size="x-large",
        ),
        bbox_to_anchor=[0.99, 0.99],
        loc="upper right",
        frameon=False,
    )
    legt.set_draggable(True)
    ax.add_artist(legt)

    if show_copyright:
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

    # Components of the solution
    lines = ["--", "-.", ":", "-"]
    linecycler = cycle(lines)

    labels = list()
    texlabels = list()

    if not (np.iterable(indices)):
        indices = (indices,)
    if not (np.iterable(offsets)):
        offsets = (offsets,)
    if len(offsets) != len(indices):
        raise AttributeError(f"lengths of {indices=} and {offsets=} do not match.")
    for i, (offset, index) in enumerate(zip(offsets, indices)):
        raw = list()
        parameters = list()
        db_idx = database_idx[index]
        db = database[db_idx]
        dbindex = index - database_off[database_idx[index]]
        if len(database) > 1:
            db_name = database_label[db_idx]
            try:
                int(db_name)
                db_name = f"DB {db_name}"
            except:
                pass
            raw.append(db_name)
            parameters.append(db_name)
        for j in range(len(db.fieldnames)):
            if db.fieldflags[j] != StarDB.Flags.parameter:
                continue
            value = db.fielddata[dbindex][j]
            unit = db.fieldunits[j]
            name = db.fieldnames[j]
            form = db.fieldformats[j]
            if unit == "-":
                unit = ""
            raw.append(f"{value:{form}} {unit}".strip())

            if name in _title_formatters:
                value = _title_formatters[name](value, form, unit)
            elif unit in _unit_formatters:
                value = _unit_formatters[unit](value, form)
            else:
                value = f"{value:{form}}"
                unit = _unit_translate.get(unit, unit)
                if unit not in (
                    "",
                    "-",
                ):
                    value = f"{value} {unit}"
            parameters.append(value)
        texlabels.append(", ".join(parameters))
        labels.append(f"{dbindex}: " + ", ".join(raw))

        # The pattern found by the algorithm
        y_a = np.log10(summed[index_t]) - logsun_full
        x_a = np.array(zlist_comb)

        x_ga = np.array(zlist_comb)

        # Move the x position for combined things
        x_star = x_star.astype("f8")
        x_a = x_a.astype("f8")
        x_ga = x_ga.astype("f8")
        for row in combine:
            if len(row) > 0:
                for x in [x_star, x_a, x_ga]:
                    # Get the first in the combined elements
                    ind = np.where(x == row[0])[0]
                    if len(ind) > 0:
                        # Make the new position the average of the first
                        # consecutives
                        consec = np.array_split(row, np.where(np.diff(row) != 1)[0] + 1)
                        x[ind[0]] = np.mean(consec[0])
        # Plot components if there are more than one
        if len(indices) > 1:
            y_ga = (
                np.log10((full_abudata[index_t, index] * offset) + 1.0e-99)
                - logsun_full
            )

            ax.plot(
                x_ga,
                y_ga,
                marker="",
                linestyle=next(linecycler),
                lw=1.3,
                color=colorsys.hsv_to_rgb(i * 360 / len(indices), 1, 1),
                # label = str(db.fielddata[index]['mass'])),
                label=texlabels[i],
            )

    # Change the label of the summed line based on how many components there are
    if len(indices) > 1:
        sumlabel = "Sum"
    else:
        sumlabel = texlabels[0]

    # Green line
    ax.plot(
        x_a,
        y_a,
        color="g",
        lw=2.0,
    )
    # Red triangles
    ax.plot(
        x_a[np.invert(lolim_index_all)],
        y_a[np.invert(lolim_index_all)],
        marker="^",
        markerfacecolor="red",
        markeredgecolor="red",
        lw=0,
    )
    # Hollow triangles
    ax.plot(
        x_a[lolim_index_all],
        y_a[lolim_index_all],
        marker="^",
        markerfacecolor="white",
        markeredgecolor="red",
        lw=0,
    )
    # Hidden line for legend purposes
    ax.plot(
        1000,
        1000,
        marker="^",
        markerfacecolor="red",
        markeredgecolor="red",
        color="g",
        lw=2.0,
        label=sumlabel,
    )

    if dist is None:
        # Calculate bounding box
        txt = fig.text(0, 0, "gh", fontsize=annosize)
        renderer = fig.canvas.get_renderer()
        bbox = txt.get_window_extent(renderer)
        txt.remove()
        dist = bbox.height

    # Plot limits
    if ylim is None:
        ylim = (np.min(y_star) - 1.0, 0.9)

    if xlim is None:
        xlim = (zlist_comb[0] - 0.99, zlist_comb[-1] + 0.99)

    # Calculate number of pixels per data
    dpi = fig.get_dpi()
    height_inch = figsize[1]
    height_pixel = height_inch * dpi
    data_scale = height_pixel / (ylim[1] - ylim[0])
    space = fontsize / data_scale
    gap_size = 3.5 * space

    anno = np.copy(y_a)

    # This is the ind corresponding to elements in the star data (the error bar points)
    for z, zor in zip(x_a, zlist_comb):
        # We use x_a instead of zlist_comb because it has the modified Z
        # for the combined elements
        ind = np.where(x_a == z)[0][0]
        if zor in range(2, zlist_comb[-1]):
            if y_a[ind] - y_a[ind - 1] >= y_a[ind + 1] - y_a[ind]:
                if z in x_star:
                    star_ind = np.where(x_star == z)[0][0]
                    if (
                        0 < (y_star[star_ind] - y_a[ind]) <= gap_size
                        or 0
                        < (y_star[star_ind] - y_a[ind] - y_star_err[1, [star_ind]])
                        <= gap_size
                    ):
                        anno[ind] = y_star[star_ind] + y_star_err[1, [star_ind]]
                loc = (0, 0.7 * dist)
            elif y_a[ind] - y_a[ind - 1] <= y_a[ind + 1] - y_a[ind]:
                if z in x_star:
                    star_ind = np.where(x_star == z)[0][0]
                    if (
                        0 < -(y_star[star_ind] - y_a[ind]) <= gap_size
                        or 0
                        < -(y_star[star_ind] - y_a[ind]) - y_star_err[0, [star_ind]]
                        <= gap_size
                    ):
                        anno[ind] = y_star[star_ind] - y_star_err[1, [star_ind]] - space
                loc = (0, -dist)

        elif z == 1 or z == zlist_comb[-1]:
            loc = (0, 0.7 * dist)

        # Make a special label for combined things
        label = None
        for row, group in enumerate(combine):
            if len(group) > 0:
                if zor == group[0]:
                    strings = [str(I(z)) for z in combine[row]]
                    label = "+".join(strings)
        if label is None:
            label = I(int(z)).element_symbol()

        ax.annotate(
            label,
            xy=(z, anno[ind]),
            xycoords="data",
            xytext=loc,
            textcoords="offset points",
            size=annosize,
            ha="center",
            va="center",
        )

    leg = ax.legend(
        bbox_to_anchor=[0.98, 0.92],
        loc="upper right",
        numpoints=1,
        prop={"size": fontsize},
        frameon=False,
    )
    leg.set_draggable(True)

    # Show detection thresholds
    for x, y in zip(x_star, y_star_det):
        if y < 20:
            continue
        ax.plot(
            x + 0.4 * np.array([-1, 1]),
            np.array([y, y]),
            ls="-",
            lw=data_size,
            color="#0000003f",
        )
    # Show correlated errors
    cov_sel = (y_star_cov > 0) & ~up_lims
    ax.errorbar(
        np.array(x_star)[cov_sel],
        y_star[cov_sel],
        yerr=y_star_uco[:, cov_sel],
        ls="None",
        marker="o",
        ms=0,
        capsize=data_size,
        color=(0.7, 0.7, 0.7),
        mfc=(0, 0, 0),
        uplims=False,
    )
    # Plot for the excluded data points
    ax.errorbar(
        np.array(x_star)[exclude_index],
        y_star[exclude_index],
        yerr=y_star_err[:, exclude_index],
        ls="None",
        marker="o",
        ms=data_size * 1.5,
        color=(0, 0, 0),
        mfc=(1, 1, 1),
        capsize=data_size,
        uplims=up_lims[exclude_index],
    )
    # Plot for the data points
    ax.errorbar(
        np.array(x_star)[~exclude_index],
        y_star[~exclude_index],
        yerr=y_star_err[:, ~exclude_index],
        ls="None",
        marker="o",
        ms=0,
        capsize=data_size,
        color=(0, 0, 0),
        mfc=(0, 0, 0),
        uplims=up_lims[~exclude_index],
    )

    # Make some dots for points that aren't uplims
    ii = np.where(~exclude_index)[0]
    ij = np.where(up_lims[ii] == 0)[0]
    ii = ii[ij]
    ax.scatter(
        np.array(x_star)[ii],
        y_star[ii],
        marker="o",
        s=10 * data_size,
        color=(0, 0, 0),
    )
    #    y_starwoBBN = y_star[2:]
    # ax.set_ybound(np.min(y_star)-1., np.max(y_starwoBBN + 1.))

    ax.set_ybound(*ylim)
    ax.set_xbound(*xlim)

    # ax.ticklabel_format(style='plain')
    ax.tick_params(axis="both", which="major", labelsize=fontsize)

    class Formatter(mpl.ticker.FuncFormatter):
        def __init__(self, *args, **kwargs):
            def formatter(*args, **kwargs):
                """
                function to format string for use with ticker
                """
                v = args[0]
                if int(v) == v:
                    v = int(v)
                return str(v)

            super().__init__(formatter)

    ax.xaxis.set_major_formatter(Formatter())
    ax.yaxis.set_major_formatter(Formatter())

    fig.tight_layout()
    fig.show()

    if savename is not None:
        if save is None:
            save = True

    if save:
        if savename is None:
            if type(indices) is list:
                savename = (
                    star.name
                    + "."
                    + ".".join([str(index) for index in indices])
                    + "."
                    + str(solution_fitness)
                    + "."
                    + save
                )
            else:
                savename = (
                    star.name
                    + "."
                    + "."
                    + str(indices)
                    + "."
                    + str(solution_fitness)
                    + "."
                    + save
                )
            savename = "../plots/" + savename + ".pdf"
        fig.savefig(savename)

    return labels, (x_a, y_a)


def fitplot(starname, generations, popsize, genesize, times, history, gen=False):
    fig = plt.figure(
        figsize=(10, 6),
        dpi=102,
        facecolor="white",
        edgecolor="white",
    )
    ax = fig.add_subplot(111)
    if gen:
        ax.set_xlabel("Generations")
        times = np.arange(generations + 1)
    else:
        ax.set_xlabel("Time (s)")
    ax.set_ylabel("Fitness (Error)")
    ax.set_title(
        "{starname:s} - Generations: {generations:d}, Population size: {popsize:d}, Gene size: {genesize:d}".format(
            starname=starname,
            generations=generations,
            popsize=popsize,
            genesize=genesize,
        )
    )
    ax.plot(times, history["average"], label="average", marker="", color="tab:blue")
    ax.plot(times, history["best"], label="best", marker="", color="tab:green")
    ax.plot(times, history["worst"], label="worst", marker="", color="tab:red")
    ax.set_yscale("log")
    leg = ax.legend(loc="best")
    leg.set_draggable(True)
    fig.tight_layout()
