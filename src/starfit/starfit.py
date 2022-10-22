"""Results objects from the various algorithms"""

from pathlib import Path
from string import capwords
from textwrap import wrap

import numpy as np

from starfit.utils import find_data

from . import DB, REF, SOLAR, starplot
from .autils import abusets
from .autils.abuset import AbuData
from .autils.isotope import ion as I
from .autils.isotope import ufunc_Z
from .autils.logged import Logged
from .autils.physconst import MSUN
from .autils.stardb import StarDB
from .fitness import solver
from .read import Star

# uniform and brief units
_unit_translate = {
    "solar masses": "Msun",
    "M_sun": "Msun",
    "He core fraction": "He core",
}

_title_translate = {
    "lower mass cut": "Mlow",
    "upper mass cut": "Mhigh",
}


class Results(Logged):
    """
    Object for running the various algorithms and a container for the results
    """

    def __init__(self):
        self.history = {"best": [], "average": [], "worst": []}
        self.initsol = None
        self.bestsol = None
        self.times = []
        self.gen = None
        self.pop_size = None

        try:
            self.silent
        except:
            self.silent = False

        self.setup_logger(silent=self.silent)

    @staticmethod
    def _normalize_mass(db, name):
        data = np.array(db.fielddata[name])
        i = db.fieldnames.tolist().index(name)
        if db.fieldunits[i] in (
            "Msun",
            "M_sun",
            "solar masses",
            "msun",
            "MSUN",
        ):
            db.fieldunits[i] = "Msun"
        elif db.fieldunits[i] in ("g",):
            db.fieldunits[i] = ("Msun",)
            data *= MSUN
        else:
            raise AttributeError(f"unknown mass unit {db.fieldunits[i]}")
        return data

    def _setup(
        self,
        filename,
        database,
        *,
        combine=None,
        z_exclude=None,
        z_min=1,
        z_max=30,
        upper_lim=None,
        z_lolim=None,
    ):
        """Prepare the data for the solvers.  Trims the databases and excludes
        elements.  Combines multiple databases.

        filename:
           star to match
        database:
           string of, path to, or object of data base.
           multiple data bases are provided as an iterable.
        """

        if isinstance(database, (str, Path, StarDB)):
            database = (database,)

        # read in data bses
        self.db = list()
        self.ejecta = list()
        self.db_n = len(database)
        self.data = list()
        for db in database:
            if not isinstance(db, StarDB):
                dbpath = find_data(DB, db)
                # Read the database
                db = StarDB(dbpath, silent=self.silent)
            self.db.append(db)

            if "ejecta" in db.fieldnames:
                ejecta = self._normalize_mass(db, "ejecta")
            elif "mass" in db.fieldnames:
                ejecta = self._normalize_mass(db, "mass")
                if "remnant" in db.fieldnames:
                    ejecta -= self._normalize_mass(db, "remnant")
            else:
                ejecta = np.full(db.fielddata.shape[0], 1.0e0, dtype=np.float64)
            self.ejecta.append(ejecta)

            data = AbuData(db.data, db.ions, molfrac=True)
            self.data.append(data)
        if self.db_n == 1:
            self.ions = self.data[0].ions.copy()
            self.data = self.data[0].data.copy()
            self.db_idx = np.full(self.data.shape[0], 0, dtype=np.int64)
            self.db_off = np.array([0], dtype=np.int64)
            self.db_num = np.array([self.data.shape[0]], dtype=np.int64)
        else:
            self.data = AbuData.join(self.data)
            self.ions = self.data.ions
            self.data = self.data.data
            self.db_idx = np.ndarray(self.data.shape[0], dtype=np.int64)
            self.db_off = np.ndarray(len(self.db), dtype=np.int64)
            self.db_num = np.ndarray(len(self.db), dtype=np.int64)
            n0 = 0
            for i, db in enumerate(self.db):
                n = db.data.shape[0]
                n1 = n0 + n
                self.db_idx[n0:n1] = i
                self.db_off[i] = n0
                self.db_num[i] = n
                n0 = n1

        if combine is None:
            combine = []
        if not np.iterable(combine):
            raise AttributeError(f"{combine=} not supported")
        if len(combine) == 0:
            combine = [[]]
        elif len(combine) == 1:
            if not np.iterable(combine[0]):
                combine = [combine]

        if z_exclude is None:
            z_exclude = []
        if upper_lim is None:
            upper_lim = []
        if z_lolim is None:
            z_lolim = []

        self.combine = combine
        self.z_min = z_min
        self.z_max = z_max
        self.upper_lim = upper_lim
        self.z_lolim = z_lolim
        self.z_exclude = z_exclude

        # Read a star
        star = Star(filename, silent=self.silent)
        self.star = star

        # Remove elements with Z > z_max and Z < z_min
        mask_zmax = np.array(
            [
                ion.Z <= z_max and ion.Z >= z_min
                for ion in star.element_abundances.element
            ]
        )
        eval_data = star.element_abundances[mask_zmax]
        # eval_data = eval_data #The other one gets chopped and changed

        # Remove upper limit elements if upper limits is not enabled
        mask_uplim = np.array([error > 0 for error in eval_data.error])
        if not upper_lim:
            eval_data = eval_data[mask_uplim]

        # List of DB ions which will get chopped
        dbions = self.ions.copy()

        # The full set of abundance data from the database
        full_abudata = self.data.copy().T
        full_ions = self.ions.copy()

        # List of every element in the DB
        list_db = np.array(self.ions)

        # List of every element in the DB less than z_max and greater than z_min
        zdb = ufunc_Z(list_db)
        list_ztrim = list_db[(zdb <= z_max) & (zdb >= z_min)]

        # Prepare the sun
        sun = abusets.SolAbu(
            name=find_data(REF, SOLAR),
            silent=self.silent,
        )

        # Transforms the uncombined element numbers for use in the sun
        index_sun = np.in1d(list_db, list_ztrim, assume_unique=True)
        sunions = dbions[index_sun]
        sun_full = sun.Y(full_ions[index_sun])
        sun_star = sun.Y(eval_data.element[:])

        # Combine elements
        list_comb = list_ztrim
        comb_abudata = full_abudata
        comb_ions = full_ions
        # Structure of list:
        #  [[6,7,8], [1,2]] will combine C+N+O and H+He
        if len(combine[0]) > 0:
            self.logger.info("Combining elements:")
            for group in combine:
                elements = [I(z) for z in group]
                for ele in elements:
                    assert (
                        ele in star.element_abundances["element"]
                    ), "Combined elements need to have an observed data point"

                self.logger.info("    " + "+".join([str(ion) for ion in elements]))
                # Get abundance of each element
                abu = [
                    10 ** eval_data.abundance[np.in1d(eval_data.element, [ele])]
                    for ele in elements
                ]

                # Replace the element with the lowest Z with the total abundance
                eval_data.abundance[
                    np.in1d(eval_data.element, [elements[0]])
                ] = np.log10(np.sum(abu))

                # Get error of each element
                err = [
                    eval_data.error[np.in1d(eval_data.element, [ele])]
                    for ele in elements
                ]

                for e in err:
                    assert e > 0, "Can only combine data points, not upper lims"

                # Replace the element with the first
                ae = [
                    (abu[i] * (10 ** np.abs(err[i]) - 1)) ** 2
                    for i in range(len(group))
                ]

                eval_data.error[np.in1d(eval_data.element, [elements[0]])] = np.log10(
                    (np.sqrt(np.sum(ae)) / np.sum(abu)) + 1
                )

                # Do the same process for the sun
                # Sum
                sunabu = [sun_full[np.in1d(sunions, [ele])] for ele in elements]
                # Write
                sun_full[np.in1d(sunions, [elements[0]])] = np.sum(sunabu)
                sun_star[np.in1d(eval_data.element[:], [elements[0]])] = np.sum(sunabu)
                # Remove
                idbse = np.invert(np.in1d(sunions, elements[1:]))
                sun_full = sun_full[idbse]
                sunions = sunions[idbse]
                iee = np.invert(np.in1d(eval_data.element[:], elements[1:]))
                sun_star = sun_star[iee]

                # Remove all the elements except the first in the group
                eval_data = eval_data[iee]

                # Do the same process for the database
                # Sum
                abudb = np.array(
                    [comb_abudata[np.in1d(dbions, [ele])][0] for ele in elements]
                )

                # Write
                comb_abudata[np.in1d(dbions, [elements[0]])][0] = np.sum(abudb, axis=0)

                # Remove
                idbe = np.invert(np.in1d(dbions, elements[1:]))
                comb_abudata = comb_abudata[idbe]
                comb_ions = comb_ions[idbe]
                dbions = dbions[idbe]

                # Adjust z list for the new combined elements
                ile = np.invert(np.in1d(list_comb, elements[1:]))
                list_comb = list_comb[ile]

        else:
            comb_abudata = full_abudata
            comb_ions = full_ions

        # Keep track of upper limits for plotting (do this after combining elements to get the indices right)
        self.uplim_index_star = eval_data.error < 0

        # Trim the DB of elements that are not in the star file
        mask_star = np.in1d(comb_ions, eval_data.element, assume_unique=True)
        trimmed_db = comb_abudata[mask_star]

        # The index of all the excluded elements
        exclude_index = np.in1d(
            [i.Z for i in eval_data.element],
            z_exclude,
            assume_unique=True,
        )

        # The index of lower limits of all elements (post-combining)
        lolim_index_all = np.in1d([i.Z for i in list_comb], z_lolim, assume_unique=True)

        # The index of lower limits of the star's elements
        lolim_index_star = np.in1d(
            [i.Z for i in eval_data.element], z_lolim, assume_unique=True
        )

        # Flip the sign of the error of the lower limits
        eval_data.error[lolim_index_star] = -np.abs(eval_data.error[lolim_index_star])

        self.logger.info(
            f"Matching {len(eval_data)} data points from Z={z_min} to Z={z_max}:"
        )

        for line in wrap(" ".join([str(ion) for ion in eval_data.element]), 60):
            self.logger.info("    " + line)

        self.logger.info(
            f"with {np.sum(eval_data.error < 0)} upper limits in the data:"
        )
        for line in wrap(
            " ".join([str(ion) for ion in eval_data.element[eval_data.error < 0]]), 60
        ):
            self.logger.info("    " + line)

        self.logger.info(f"and {np.sum(lolim_index_star)} lower limits in the models:")
        for line in wrap(
            " ".join([str(ion) for ion in eval_data.element[lolim_index_star]]), 60
        ):
            self.logger.info("    " + line)

        self.eval_data = eval_data
        self.trimmed_db = trimmed_db
        self.full_abudata = full_abudata
        self.full_ions = full_ions
        self.list_db = list_db
        self.list_ztrim = list_ztrim
        self.list_comb = list_comb
        self.exclude_index = exclude_index
        self.lolim_index_all = lolim_index_all
        self.lolim_index_star = lolim_index_star
        self.sun = sun
        self.sun_full = sun_full
        self.sun_star = sun_star
        self.fit_size = self.trimmed_db.shape[0]
        self.db_size = self.trimmed_db.shape[1]

        del self.data
        del self.ions

    def run(
        self,
        stars,
        offsets=None,
        fixed=False,
    ):
        """Solve for the specified stars"""
        stars = np.array(stars)
        sol = np.ndarray(stars.shape, dtype=[("index", "int"), ("offset", "f8")])
        for i in range(stars.shape[0]):
            sol[i, :]["index"] = stars[i, :]
            if offsets is None:
                sol[i, :]["offset"] = 1.0e-3 / stars.shape[1]
                ls = True
            else:
                sol[i, :]["offset"] = offsets[i]
                ls = False

        fitness = _fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            sol,
            self.ejecta,
            fixed_offsets=fixed,
            cdf=self.cdf,
            ls=ls,
        )
        return sol, fitness

    def plot(self, index=0, **kwargs):
        """Call plotting routines to plot the best fit."""

        bestsol = self.sorted_stars[index]
        self.labels, self.plotdata = starplot.abuplot(
            indices=bestsol["index"].tolist(),
            offsets=bestsol["offset"].tolist(),
            star=self.star,
            database=self.db,
            database_idx=self.db_idx,
            database_off=self.db_off,
            full_abudata=self.full_abudata,
            eval_data=self.eval_data,
            list_db=self.list_db,
            list_comb=self.list_comb,
            sun_full=self.sun_full,
            sun_star=self.sun_star,
            combine=self.combine,
            exclude_index=self.exclude_index,
            uplim_index_star=self.uplim_index_star,
            lolim_index_all=self.lolim_index_all,
            solution_fitness=self.sorted_fitness[0],
            **kwargs,
        )

    def plot_fitness(self):
        # Fitness over time plot
        raise NotImplementedError()

    def text_result(self, n=20, *, n0=0, format="unicode", wide=12):
        """Print data of best fit."""
        text = []
        if format == "html":
            base_title = ["&#x1D6D8;&sup2;", "Dilution", "Index"]
        elif format == "unicode":
            base_title = ["\u03C7\u00B2", "Dilution", "Index"]
        else:
            base_title = ["chi**2", "Dilution", "Index"]
        base_units = ["", "(log)", ""]

        if self.db_n > 1:
            base_title[2:2] = ["DB"]
            base_units[2:2] = [""]
        empty_title = [""] * len(base_title)

        n1 = min(n0 + n, len(self.sorted_stars))

        dbidx = np.array(self.db_idx[self.sorted_stars[n0:n1]["index"]])
        dbx = np.unique(dbidx.flat)

        # wide format
        x = [
            [(fn, fu) for fn, fu in zip(db.fieldnames, db.fieldunits)]
            for db in [self.db[i] for i in dbx]
        ]
        ii = np.full((len(x), np.max([len(_) for _ in x])), -1, dtype=np.int64)
        fields = list()
        nfield = 0
        for j, y in enumerate(x):
            for jj, yy in enumerate(y):
                if yy in fields:
                    ii[j, jj] = fields.index(yy)
                else:
                    ii[j, jj] = nfield
                    fields.append(yy)
                    nfield += 1
        fieldmap = np.full((self.db_n, ii.shape[1]), -1, dtype=np.int64)
        for j, i in enumerate(dbx):
            fieldmap[i] = ii[j]
        fields = np.array(fields)

        if isinstance(wide, int):
            wide = len(base_title) + nfield <= wide

        def _head(base_title, base_units, title, units):
            text.append(
                base_title
                + [
                    _title_translate[word]
                    if word in _title_translate
                    else capwords(word)
                    for word in title
                ]
            )
            text.append(
                base_units
                + [
                    f"({_unit_translate.get(word, word)})"
                    if word
                    not in (
                        "",
                        "-",
                    )
                    else ""
                    for word in units
                ]
            )

        def _wide_head():
            _head(base_title, base_units, fields[:, 0], fields[:, 1])

        def _short_head(db, full=False):
            if full:
                title = base_title
                units = base_units
            else:
                title = units = empty_title
            _head(title, units, db.fieldnames, db.fieldunits)

        db_idx0 = dbidx[0, 0]
        if wide:
            _wide_head()
        else:
            _short_head(self.db[db_idx0], True)

        for i in range(n0, n1):
            for j in range(self.sol_size):
                index, offset = self.sorted_stars[i, j]
                db_idx = self.db_idx[index]
                db = self.db[db_idx]
                if db_idx != db_idx0 and not wide:
                    _short_head(db, j == 0 and self.sol_size > 1)
                    db_idx0 = db_idx
                dbindex = index - self.db_off[db_idx]
                data = db.fielddata[dbindex]
                line = list()
                if j == 0:
                    line.append(f"{self.sorted_fitness[i]:3.2f}")
                else:
                    line.append("")
                line.append(f"{np.log10(offset):7.2f}")
                if self.db_n > 1:
                    line.append(f"{db_idx + 1:>d}")
                line.append(f"{dbindex:6d}")
                if wide:
                    for k in range(nfield):
                        if k in fieldmap[db_idx]:
                            l = np.where(fieldmap[db_idx] == k)[0][0]
                            line.append(f"{data[l]:{db.fieldformats[l]}}")
                        else:
                            line.append("")
                else:
                    line.extend([f"{x:{y}}" for x, y in zip(data, db.fieldformats)])
                text.append(line)
            if self.sol_size > 1:
                text.append("")

        return text

    @staticmethod
    def textpad(s, n):
        nspace = n - len(s)
        return " " * nspace + s

    def print(self, *args, **kwargs):
        full = kwargs.pop("full", False)
        print(self.format(*args, **kwargs))
        if full:
            self.print_comments()
        else:
            print(self.format_db())

    def format_db(self, ind=0):
        pad = " " * ind
        string = [pad + "DB  Name"]
        for i, db in enumerate(self.db):
            string.append(pad + f"{i+1:>2d}  {db.name}")
        string = "\n".join(string)
        return string

    def print_db(self, ind=0):
        print(self.format_db(ind))

    def format_comments(self, npad=72):
        string = list()
        for i, db in enumerate(self.db):
            if i == 0:
                string.append("=" * npad)
            if self.db_n > 1:
                string.append(f"{i+1:>2d}:  {db.name}")
            else:
                string.append(f"{db.name}")
            string.append("-" * npad)
            string.extend(db.comments.tolist())
            string.append("=" * npad)
        string = "\n".join(string)
        return string

    def print_comments(self, ind=0):
        print(self.format_comments())

    def format(self, *args, **kwargs):
        text = self.text_result(*args, **kwargs)
        lengths = [[len(word) for word in line] for line in text if len(line) > 0]
        lengths = np.asarray(lengths).max(axis=0) + 1

        string = ""
        for line in text:
            for word, length in zip(line, lengths):
                string += self.textpad(word, length)
            string += "\n"
        return string

    def __str__(self):
        return self.format(n=10)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.star.name})"


def _fitness(
    trimmed_db,
    eval_data,
    z_exclude_index,
    sol,
    ejecta=[],
    fixed_offsets=False,
    cdf=0,
    ls=False,
):
    """
    Evaluate the fitness of a set of solutions.
    If abundance is directly given in the case of smart GA, use that.
    """

    if fixed_offsets:
        offset = ejecta[sol["index"]]
        ls = False
    else:
        offset = sol["offset"]

    error = np.copy(eval_data.error)
    error[z_exclude_index] = 1.0e99

    abu = np.transpose(trimmed_db[:, sol["index"]], (1, 2, 0))

    fitness, offsets = solver.fitness(
        eval_data.abundance,
        error,
        abu,
        offset,
        cdf=cdf,
        ls=ls,
    )

    sol["offset"] = offsets
    fitness /= eval_data.error.shape[0] - np.sum(z_exclude_index) - 1
    return fitness
