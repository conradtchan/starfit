"""Results objects from the various algorithms"""

from pathlib import Path
from string import capwords
from textwrap import wrap

import numpy as np

from . import DB, REF, SOLAR
from .autils import abusets
from .autils.abuset import AbuData
from .autils.isotope import ion as I
from .autils.isotope import ufunc_Z
from .autils.logged import Logged
from .autils.physconst import MSUN
from .autils.stardb import StarDB
from .autils.utils import index1d, is_iterable
from .fit import get_fitness
from .star import Star
from .starplot import abuplot
from .utils import find_all, find_data

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


DB_LABEL_SIZE = 8


class StarFit(Logged):
    """
    Object for running the various algorithms and a container for the results
    """

    def __init__(self, *args, silent=False, **kwargs):
        self.history = {"best": [], "average": [], "worst": []}
        self.initsol = None
        self.bestsol = None
        self.times = []
        self.gen = None
        self.pop_size = None
        self.silent = silent

        self.setup_logger(silent=self.silent)

        self._setup(*args, **kwargs)

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
        db,
        *,
        combine=None,
        z_exclude=None,
        z_min=1,
        z_max=999,
        upper_lim=True,
        z_lolim=None,
        y_floor=1.0e-99,
        db_label=None,
        cdf=True,
        cov=False,
        debug=False,
        show=False,
    ):
        """Prepare the data for the solvers.  Trims the databases and excludes
        elements.  Combines multiple databases.

        filename:
           star to match
        database:
           string of, path to, or object of data base.
           multiple data bases are provided as an iterable.
        """
        self.debug = debug

        if isinstance(db, (str, Path, StarDB)):
            db = (db,)
        database = list()
        for d in db:
            if d in ("*", ..., "...", "**"):
                database.extend(
                    find_all(
                        DB,
                        "*.stardb.*",
                        complete=d
                        in (
                            Ellipsis,
                            "...",
                            "**",
                        ),
                    )
                )
            elif str(d).find("*") >= 0 or str(d).find("*") >= 0:
                database.extend(find_all(DB, d))
            else:
                database.append(d)
        assert len(database) > 0, f"require valid database, {db=}"

        # read in databases
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
            self.data = AbuData.join(self.data, silent=True)
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

        if db_label is None:
            self.db_lab = list()
            for i, d in enumerate(self.db):
                if hasattr(d, "label"):
                    self.db_lab.append(d.label)
                else:
                    self.db_lab.append(f"{i:d}")
        else:
            self.db_lab = db_label.copy()
        assert (
            is_iterable(self.db_lab) and len(self.db_lab) == self.db_n
        ), "number of labels ({len(self.db_lab)}) does not match number of databases ({self.db_n})."
        for i, l in enumerate(self.db_lab):
            assert isinstance(l, str), f"require labels to be strings: '{l}'"
            if len(l) > DB_LABEL_SIZE:
                self.logger.warn(
                    f"truncating label to {DB_LABEL_SIZE:d} chracters: '{l}' -> '{l[:4]}'"
                )
                self.db_lab[i] = l[:DB_LABEL_SIZE]

        if combine is None:
            combine = []
        if not is_iterable(combine):
            raise AttributeError(f"{combine=} not supported")
        if len(combine) == 0:
            combine = [[]]
        elif not is_iterable(combine[0]):
            combine = [combine]
        for i, group in enumerate(combine):
            for j, element in enumerate(group):
                if isinstance(element, str):
                    group[j] = I(element, element=True).Z
            combine[i] = sorted(group)
        if z_exclude is None:
            z_exclude = []
        for i, element in enumerate(z_exclude):
            if isinstance(element, str):
                z_exclude[i] = I(element, element=True).Z
        if z_lolim is None:
            z_lolim = []
        for i, element in enumerate(z_lolim):
            if isinstance(element, str):
                z_lolim[i] = I(element, element=True).Z
        if z_min is None:
            z_min = 1
        if isinstance(z_min, str):
            z_min = I(z_min, element=True).Z
        if z_max is None:
            z_max = 92
        if isinstance(z_max, str):
            z_max = I(z_max, element=True).Z

        self.z_min = z_min
        self.z_max = z_max
        self.combine = combine
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
        if cov is False:
            data_type = star.data_type.copy()
            assert data_type[-1][0] == "covariance"
            data_type[-1] = data_type[-1][:2] + (0,)

            eval_data_new = np.recarray(eval_data.shape[0], dtype=data_type)
            eval_data_new.element = eval_data.element
            eval_data_new.abundance = eval_data.abundance
            eval_data_new.detection = eval_data.detection
            eval_data_new.error = np.sign(eval_data.error) * np.sqrt(
                eval_data.error**2 + np.sum(eval_data.covariance**2, axis=1)
            )
            eval_data = eval_data_new

        # Remove upper limit elements if upper limits is not enabled
        mask_uplim = np.array([error > 0 for error in eval_data.error])
        if not upper_lim:
            eval_data = eval_data[mask_uplim]

        # List of DB ions which will get chopped
        dbions = self.ions.copy()

        # The full set of abundance data from the database
        full_abudata = self.data.copy().T
        full_ions = self.ions.copy()

        # set floor value for data
        if y_floor is not None:
            full_abudata = np.maximum(full_abudata, y_floor)

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
                elements = [I(z, element=True) for z in group]
                for element in elements:
                    assert (
                        element in star.element_abundances["element"]
                    ), "Combined elements need to have an observed data point"

                self.logger.info("    " + "+".join([str(ion) for ion in elements]))
                # Get abundance of each element
                ii = index1d(elements, eval_data.element)
                abu = 10 ** eval_data.abundance[ii]

                # Replace the element with the lowest Z with the total abundance
                i0 = np.searchsorted(eval_data.element, elements[0])
                eval_data.abundance[i0] = np.log10(np.sum(abu))

                # Get error of each element
                error = eval_data.error[ii]

                assert np.all(error > 0), "Can only combine data points, not upper lims"

                # Replace the first element with the combined value
                ae = abu * (10**error - 1) ** 2
                eval_data.error[i0] = np.log10(np.sqrt(np.sum(ae) / np.sum(abu)) + 1)

                # Do the same for detection thresholds
                det = eval_data.detection[ii]
                ad = abu[:, np.newaxis] * (10**det - 1)
                eval_data.detection[i0] = np.log10(np.sum(ad) / np.sum(abu) + 1)

                # Do the same for covariances
                cov = eval_data.covariance[ii]
                ac = abu[:, np.newaxis] * (10**cov - 1)
                eval_data.covariance[i0] = np.log10(
                    np.sum(ac, axis=0) / np.sum(abu) + 1
                )

                # Do the same process for the sun
                # Sum
                sunabu = np.sum(sun_full[index1d(elements, sunions)])
                # Write
                sun_full[np.searchsorted(sunions, elements[0])] = sunabu
                sun_star[i0] = sunabu

                # Remove
                idbse = np.invert(np.in1d(sunions, elements[1:]))
                sun_full = sun_full[idbse]
                sunions = sunions[idbse]
                iee = np.invert(np.in1d(eval_data.element, elements[1:]))
                sun_star = sun_star[iee]

                # Remove all the elements except the first in the group
                eval_data = eval_data[iee]

                # Do the same process for the database
                # Sum
                abudb = np.sum(comb_abudata[index1d(elements, dbions)], axis=0)

                # Write
                comb_abudata[np.searchsorted(sunions, elements[0])] = abudb

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

        n = np.count_nonzero(eval_data.error < 0)
        if n > 0:
            self.logger.info(f"including {n} upper limits in the data:")
            for line in wrap(
                " ".join([str(ion) for ion in eval_data.element[eval_data.error < 0]]),
                60,
            ):
                self.logger.info("    " + line)

        n = np.count_nonzero(lolim_index_star)
        if n > 0:
            self.logger.info(f"including {n} lower limits in the models:")
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
        self.cdf = cdf

        del self.data
        del self.ions

        if show:
            print()
            self.print_db(ind=4)
        self.show = show

    def run(
        self,
        stars,
        offsets=None,
        fixed_offsets=False,
        optimize=True,
    ):
        """Solve for the specified stars"""
        stars = np.array(stars)
        sol = np.ndarray(stars.shape, dtype=[("index", "int"), ("offset", "f8")])
        for i in range(stars.shape[0]):
            sol[i, :]["index"] = stars[i, :]
            if offsets is None:
                assert optimize is True
                sol[i, :]["offset"] = 1.0e-3 / stars.shape[1]
                ls = True
            else:
                sol[i, :]["offset"] = offsets[i]
                if optimize is True:
                    ls = False
                else:
                    ls = None

        fitness = get_fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            sol,
            self.ejecta,
            fixed_offsets=fixed_offsets,
            cdf=self.cdf,
            ls=ls,
        )
        return sol, fitness

    def plot(self, index=0, **kwargs):
        """Call plotting routines to plot the best fit."""

        bestsol = self.sorted_stars[index]
        self.labels, self.plotdata = abuplot(
            indices=bestsol["index"].tolist(),
            offsets=bestsol["offset"].tolist(),
            star=self.star,
            database=self.db,
            database_idx=self.db_idx,
            database_off=self.db_off,
            database_label=self.db_lab,
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

    def text_result(self, n=20, *, n0=0, format="unicode", wide=12, _return_dbx=False):
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
        nfield_short = 0
        for j, y in enumerate(x):
            for jj, yy in enumerate(y):
                if yy in fields:
                    ii[j, jj] = fields.index(yy)
                else:
                    ii[j, jj] = nfield
                    fields.append(yy)
                    nfield += 1
            nfield_short = max(nfield_short, len(y))
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
            pad = [""] * (nfield_short - db.nfield)
            fieldnames = list(db.fieldnames) + pad
            fieldunits = list(db.fieldunits) + pad
            _head(title, units, fieldnames, fieldunits)

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
                    line.append(f"{self.db_lab[db_idx]}")
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
                    line.extend([""] * (nfield_short - db.nfield))
                text.append(line)
            if self.sol_size > 1:
                text.append("")

        if _return_dbx:
            return text, dbx
        return text

    @staticmethod
    def textpad(s, n):
        nspace = n - len(s)
        return " " * nspace + s

    def print(self, *args, **kwargs):
        full = kwargs.pop("full", False)
        kwargs["_return_dbx"] = True
        text, dbx = self.format(*args, **kwargs)
        print(text)
        if self.db_n == 1:
            return
        if full is True:
            self.print_comments(dbx=dbx)
        elif full is False:
            self.print_db(dbx=dbx)

    def format_db(self, ind=0, pad="", dbx=None):
        if pad is None:
            pad = ""
        pad = " " * ind + pad
        if dbx is None:
            dbx = range(self.db_n)
        db_len = 2
        num = True
        for i in dbx:
            lab = self.db_lab[i]
            db_len = max(db_len, len(lab))
            try:
                num &= int(lab) == i
            except:
                pass
        if num:
            string = [pad + f"{'DB':>{db_len}}  Name"]
        else:
            string = [pad + f"NR  {'DB':>{db_len}}  Name"]
        for i in dbx:
            db = self.db[i]
            if num:
                string.append(pad + f"{self.db_lab[i]:>{db_len}}  {db.name}")
            else:
                string.append(pad + f"{i:>2d}  {self.db_lab[i]:>{db_len}}  {db.name}")
        string = "\n".join(string)
        return string

    def print_db(self, *args, **kwargs):
        logger = kwargs.pop("logger", False)
        text = self.format_db(*args, **kwargs)
        if logger:
            for l in text.splitlines():
                self.logger.info(l)
            return
        print(text)

    def format_comments(self, npad=72, dbx=None):
        string = list()
        if dbx is None:
            dbx = range(self.db_n)
        for i in dbx:
            db = self.db[i]
            if i == 0:
                string.append("=" * npad)
            if self.db_n > 1:
                string.append(f"{self.db_lab[i]}:  {db.name}")
            else:
                string.append(f"{db.name}")
            string.append("-" * npad)
            string.extend(db.comments.tolist())
            string.append("=" * npad)
        string = "\n".join(string)
        return string

    def print_comments(self, *args, **kwargs):
        print(self.format_comments(*args, **kwargs))

    def format(self, *args, **kwargs):
        _return_dbx = kwargs.get("_return_dbx", False)
        result = self.text_result(*args, **kwargs)
        if _return_dbx:
            text, dbx = result
        else:
            text = result
        lengths = [[len(word) for word in line] for line in text if len(line) > 0]
        lengths = np.asarray(lengths).max(axis=0) + 1

        string = ""
        for line in text:
            for word, length in zip(line, lengths):
                string += self.textpad(word, length)
            string += "\n"
        if _return_dbx:
            return string, dbx
        return string

    def info(self, i=0, **kwargs):
        kwargs["n0"] = i
        kwargs["n"] = 1
        kwargs.setdefault("wide", 20)
        self.print(**kwargs)

    def __str__(self):
        return self.format(n=10)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.star.name})"

    def _setup_group(self, group):
        if group is True:
            group = [[i] for i in range(self.db_n)]
        elif group is False:
            group = [[i for i in range(self.db_n)]]
        elif is_iterable(group):
            if np.all([isinstance(i, int) for i in group]):
                gnew = list()
                ndb = 0
                for nmembers in group:
                    assert nmembers > 0, "invalid group size {nmembers}"
                    gnew.append([idb + ndb for idb in range(nmembers)])
                    ndb += nmembers
                assert (
                    ndb == self.db_n
                ), f"invalid total number of databses in group {ndb=}, vs. db_n={self.db_n}"
                group = gnew
            else:
                all_db = list()
                for members in group:
                    assert is_iterable(members)
                    for idb in members:
                        assert isinstance(idb, int)
                        assert idb >= 0
                        assert idb < self.db_n
                        assert idb not in all_db
                        all_db.append(idb)
                assert len(all_db) == self.db_n
                assert len(set(all_db)) == self.db_n
        else:
            raise AttributeError(f"[{self.__class__.__name__}] unknown {group=}")

        group_n = len(group)
        group_num = np.array([np.sum(self.db_num[g]) for g in group])

        group_index = np.ndarray(self.db_size, dtype=np.int64)
        group_idx = np.ndarray(self.db_n, dtype=np.int64)
        i0 = 0
        for i, members in enumerate(group):
            for idb in members:
                n = self.db_num[idb]
                i1 = i0 + n
                group_index[i0:i1] = np.arange(n) + self.db_off[idb]
                i0 = i1
                group_idx[idb] = i

        group_off = np.cumsum(group_num, dtype=np.int64)
        group_off[1:] = group_off[0:-1]
        group_off[0] = 0

        self.group = group
        self.group_n = group_n
        self.group_num = group_num
        self.group_off = group_off
        self.group_idx = group_idx
        self.group_index = group_index
