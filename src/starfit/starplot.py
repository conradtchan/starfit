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

from .autils.isotope import Ion
from .autils.keputils import mass_string


def abuplot(
    indices=None,
    offsets=None,
    star=None,
    db=None,
    full_abudata=None,
    eval_data=None,
    list_db=None,
    list_comb=None,
    sun_full=None,
    sun_star=None,
    combine=None,
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

    x_star = np.array([ion.Z for ion in eval_data.element])
    y_star_err = np.tile(eval_data.error, (2, 1))
    # Determine upper limit points
    up_lims = uplim_index_star
    y_star_err[:, up_lims] *= np.array([-1, 0])[:, np.newaxis]

    ax.text(
        0.96,
        0.95,
        star.name,
        size="x-large",
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
    )

    if show_copyright:
        ax.text(
            0.01,
            0.01,
            r"$\copyright$ www.StarFit.org",
            size="x-small",
            horizontalalignment="left",
            verticalalignment="bottom",
            transform=ax.transAxes,
        )

    # Components of the solution
    lines = ["--", "-.", ":", "-"]
    linecycler = cycle(lines)

    labels = []
    texlabels = []

    if type(indices) is list:
        n_sol = len(indices)
    else:
        n_sol = 1
    for i in range(0, n_sol):
        try:
            index = indices[i]
            offset = offsets[i]
        except:
            index = indices
            offset = offsets

        fielddata = np.copy(db.fielddata)
        field_names = list(fielddata.dtype.names)
        indexstring = f"{index}: "
        labelfields = []
        for col, name in enumerate(field_names):
            if name == "energy":
                numformat = "{:4.2f}"
            elif name == "mixing":
                value = fielddata[index][col]
                if value != 0:
                    fielddata[index][col] = np.log10(value)
                    name = "log(mixing)"
                    nomix = False
                else:
                    nomix = True
                numformat = "{:6.4f}"
            elif name == "offset":
                numformat = "{:6.5f}"
            else:
                numformat = "{}"
            labelfields += [("{} = " + numformat).format(name, fielddata[index][col])]

        labels += [indexstring + ", ".join(labelfields) + f", dilution = {offset:4.2f}"]

        if not nomix:
            texlabels += [
                r"${} \mathrm{{M}}_\odot$, ${} \mathrm{{B}}$, $\log(f_\mathrm{{mix}})={}$".format(
                    mass_string(fielddata[index]["mass"]),
                    round(fielddata[index]["energy"], 2),
                    round(fielddata[index]["mixing"], 4),
                )
            ]
        else:
            texlabels += [
                r"${} \mathrm{{M}}_\odot$, ${} \mathrm{{B}}$, $\mathrm{{no\ mixing}}$".format(
                    mass_string(fielddata[index]["mass"]),
                    round(fielddata[index]["energy"], 2),
                )
            ]

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
        if n_sol != 1:
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
        xlim = (0.1, zlist_comb[-1] + 1)

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
                    strings = [str(Ion(z)) for z in combine[row]]
                    label = "+".join(strings)
        if label is None:
            label = Ion(int(z)).element_symbol()

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

    ax.legend(
        bbox_to_anchor=[0.98, 0.92],
        loc="upper right",
        numpoints=1,
        prop={"size": fontsize},
    ).draw_frame(False)

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

    # What is zlabels?
    return labels, (x_a, y_a)


def fitplot(sol_type, starname, generations, popsize, genesize, times, history):
    fig = plt.figure(
        figsize=(10, 6),
        dpi=102,
        facecolor="white",
        edgecolor="white",
    )
    if sol_type in ["Single", "Double"]:
        raise TypeError("Plot not available for this type of solution")
    else:
        ax_fit = fig.add_subplot(111)
        ax_fit.set_xlabel("Time (s)")
        ax_fit.set_ylabel("Fitness (Error)")
        ax_fit.set_title(
            "{starname:s} - Generations: {generations:d}, Population size: {popsize:d}, Gene size: {genesize:d}".format(
                starname=starname,
                generations=generations,
                popsize=popsize,
                genesize=genesize,
            )
        )
        ax_fit.plot(
            times, history["average"], label="average", marker="", color="tab:blue"
        )
        ax_fit.plot(times, history["best"], label="best", marker="", color="tab:green")
        ax_fit.plot(times, history["worst"], label="worst", marker="", color="tab:red")
        ax_fit.set_yscale("log")
        ax_fit.legend()
        fig.tight_layout()

    return
