"""Results objects from the various algorithms"""

import colorsys
from io import BytesIO
from itertools import cycle
from pathlib import Path
from string import capwords
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma

from . import DB
from .autils.abuset import AbuData
from .autils.isotope import Ion
from .autils.isotope import ion as I
from .autils.isotope import ufunc_idx, ufunc_Z
from .autils.logged import Logged
from .autils.physconst import MSUN
from .autils.stardb import StarDB
from .autils.utils import index1d, is_iterable
from .fit import get_fitness
from .fitness.solver import get_complete_inverse, get_complete_matrix
from .star import LOW, Star
from .starplot import (
    IntFormatter,
    _plot_title_formatters,
    _plot_unit_formatters,
    _plot_unit_translate,
    leg_copyright,
    leg_starname,
)
from .utils import find_all, find_data

# uniform and brief units
_unit_translate = {
    "solar masses": "Msun",
    "M_sun": "Msun",
    "He core fraction": "He core",
}

_unit_translate_html = {
    "Msun": "M<sub>&#x2609;</sub>",
}

_title_translate = {
    "lower mass cut": "Mlow",
    "upper mass cut": "Mhigh",
}

_title_translate_html = {
    "Mlow": "M<sub>low</sub>",
    "Mhigh": "M<sub>high</sub>",
}


DB_LABEL_SIZE = 8


class ConstraintsParseError(Exception):
    pass


class StarFit(Logged):
    r"""
    Base class for running the various fitting algorithms.

    Args:
        filename (str): Filename of star. Can be absolute or relative path. The
            files will also be searched for in the distribution files and in the
            search path specified by environment variable ``STARFIT_DATA``
            in subdirectory ``stars``.
        db (str or :class:`pathlib.Path`, optional): database file or tuple of
            data base files.  String or ``Path`` object.  Can be absolute or
            relative path.  Files will also be searched in the distribution
            files and in the search path specified by environment variable
            ``STARFIT_DATA`` in subdirectory ``db``.  You may also use the
            wildcard (``"*"``) in the data base name.  The code will then try to
            resolve all matching data bases in the first source directory that
            contains any matching file. The plain ``*`` argument will include
            all data bases in the first source that contains any data base; the
            matching is done against the pattern ``*.stardb.*``.  The
            ``Ellipis`` (``...`` Python object, not in quotation marks) will do
            the same as the plain ``*`` argument, but will continue searching
            through all data souces.  This allows for an easy way to search
            across all model data bases available. combine (list, optional): A
            list of lists of element charge numbers to treat as combined
            abundances (e.g. combine the CNO elements using ``[[6,7,8]]``).
        z_min (int, optional): Lowest element charge number to fit.
        z_max (int, optional): Highest element charge number to fit.
        z_exclude (list, optional): Element charge numbers to exclude from fit.
        lim_exclude (bool, optional): Treat ``z_min`` and ``z_max`` limits as
            exclusions (default: ``True``).  Otherwise databases are "trimmed"
            to save memory and data cannot be plotted in interactive mode.
        z_lolim (list, optional): Elements that are *model* lower limits
            (effectively the same as *observational* upper limits).
        upper_lim (bool, optional): Include observational upper limits in data
            fitting.
        cdf (bool, optional): Use the uncertainty of upper limits to calculate a
            cumulative distribution function when calculating error contribution
            (otherwise treat the upper limit as a simple one-sided
            :math:`{\chi^2}` error).
        det (bool, optional): Use the detection limits when calculating error
            contribution (experimental).
        cov (bool, optional): Use the error covariances when calculating error
            contribution (experimental).
        dst (bool, optional): Use statistical error only for detection treshold
            (default: True; experimental).
        limit_solver (bool, optional): Solver/search will only allow solutions
            for each star that contribute no more than 100%.
        limit_solution (bool, optional): Solver/search will only allow solutions
            where the total adds up to no more than 100% contributions from all
            stars.  Results from the search are renormalised accordingly.
        y_floor (bool, optional): Floor value for abundances to assume in models
            (default: ``1e.0e-99``).  This is useful for elements not produced
            in a model, otherwise :math:`{\chi^2}`; or :math:`{-\infty}`; may result.
        db_label (list, optional): A list of labels for the data bases to be
            used in plots and tables.  Will only be shown if there is more than
            one database specified.  If present, needs to match the number of
            databases specified.  If not present, databases will be numbered
            starting with ``0``, unless the ``StarDB`` has a ``label`` field
            that will be used instead.  The maximum label length currently
            allowed is ``8``.
        show (bool, optional): Show list of loaded databases with label/number
        and name, then quit.
        constraints (str, optional): String with list of conditions separated by
            comma (acts as "and").  Conditions for specific databases can be
            prefixed with a number (zero-based) if the index of the database in
            the list, followed by a colon(``:``).  Entries for different
            databases are separated by semicolon (``;``).  The fieldname has to
            be given first, then the operator, and finally the comparison value.
            Allowed operators are ``<``, ``<=``, ``==``, ``>=``, ``>``, and
            ``!=``.
        constraints_error (str, optional): one of ``warn`` (default), ``raise``,
            ``ignore``.  How StarDB deal with errors in ``constraints``.

    """

    def __init__(self, *args, silent=False, **kwargs):
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
        lim_exclude=True,
        upper_lim=True,
        z_lolim=None,
        y_floor=1.0e-99,
        db_label=None,
        cdf=True,
        cov=False,
        det=False,
        dst=False,
        debug=False,
        show=False,
        limit_solution=None,
        limit_solver=None,
        constraints=None,
        constraints_error="warn",
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

        # Read a star
        star = Star(filename, silent=self.silent)
        self.star = star

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

        constraintx = {i: list() for i in range(len(database))}
        self.constraints_ok = True
        try:
            if constraints is not None:
                for c in constraints.split(";"):
                    cx = [x.strip() for x in c.strip().split(":")]
                    if len(cx) > 2 or len(cx) == 0:
                        raise ConstraintsParseError(constraints)
                    if len(cx) == 2:
                        try:
                            idb = int(cx[0])
                            assert idb >= 0 and idb < len(database)
                        except:
                            raise ConstraintsParseError(constraints)
                        constraintx[idb].append(cx[1])
                    else:
                        for v in constraintx.values():
                            v.append(cx[0])
        except:
            self.constraints_ok = False
        for k, v in constraintx.items():
            if len(v) == 0:
                constraintx[k] = None
            else:
                constraintx[k] = ", ".join(v)

        # read in databases
        self.db = list()
        self.ejecta = list()
        self.db_n = len(database)
        self.data = list()
        for db, c in zip(database, constraintx.values()):
            if not isinstance(db, StarDB):
                dbpath = find_data(DB, db)
                # Read the database
                db = StarDB(dbpath, silent=self.silent)
            if c is not None:
                try:
                    db = db.select(c, error=constraints_error)
                except:
                    self.constraints_ok = False
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
        if not self.constraints_ok:
            self.logger.info("*** Constraints failed. ***")
        if self.db_n == 1:
            self.data = self.data[0]
            self.ejecta = self.ejecta[0]
        else:
            self.data = AbuData.join(self.data, silent=self.silent)
            self.ejecta = np.concatenate(self.ejecta)
        if not (
            set(ufunc_idx(self.star.element_abundances.element))
            < set(ufunc_idx(self.data.ions))
        ):
            star_ions = AbuData(
                np.ndarray((0, self.star.n_elements)),
                self.star.element_abundances.element,
                molfrac=True,
            )
            self.data = AbuData.join([self.data, star_ions], silent=self.silent)
        self.ions = self.data.ions
        self.data = self.data.data
        if self.db_n == 1:
            self.db_idx = np.full(self.data.shape[0], 0, dtype=np.int64)
            self.db_off = np.array([0], dtype=np.int64)
            self.db_num = np.array([self.data.shape[0]], dtype=np.int64)
        else:
            self.db_idx = np.ndarray(self.data.shape[0], dtype=np.int64)
            self.db_off = np.ndarray(self.db_n, dtype=np.int64)
            self.db_num = np.ndarray(self.db_n, dtype=np.int64)
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
                if hasattr(d, "label") and len(d.label) > 0:
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
            if isinstance(element, Ion):
                z_exclude[i] = element.Z
        if z_lolim is None:
            z_lolim = []
        for i, element in enumerate(z_lolim):
            if isinstance(element, str):
                z_lolim[i] = I(element, element=True).Z
            if isinstance(element, Ion):
                z_lolim[i] = element.Z
        if z_min is None:
            z_min = 1
        if isinstance(z_min, str):
            z_min = I(z_min, element=True).Z
        if isinstance(z_min, Ion):
            z_min = z_min.Z
        z_min = max(z_min, 1)
        if z_max is None:
            z_max = 92
        if isinstance(z_max, str):
            z_max = I(z_max, element=True).Z
        if isinstance(z_max, Ion):
            z_max = z_max.Z
        if z_max == 0:
            z_max = 92

        self.z_min = z_min
        self.z_max = z_max
        self.combine = combine
        self.upper_lim = upper_lim
        self.z_lolim = z_lolim
        self.z_exclude = z_exclude
        self.limit_solution = limit_solution
        self.limit_solver = limit_solver
        self.dst = dst

        if lim_exclude:
            z_min_ = z_min
            z_min = min(self.ions[0].Z, self.star.element_abundances.element[0].Z)
            z_exclude = list(range(z_min, z_min_)) + z_exclude
            z_max_ = z_max
            z_max = max(self.ions[-1].Z, self.star.element_abundances.element[-1].Z)
            z_exclude = z_exclude + list(range(z_max_ + 1, z_max + 1))

        # Remove elements with Z > z_max and Z < z_min
        mask_zmax = np.array(
            [
                ion.Z <= z_max and ion.Z >= z_min
                for ion in star.element_abundances.element
            ]
        )
        eval_data = star.element_abundances[mask_zmax]

        # check whether to use covariances
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

        # check whether to use detection thresholds
        if det is False:
            eval_data.detection = LOW[star.data_format]

        # Remove upper limit elements if upper limits is not enabled
        if not self.upper_lim:
            eval_data = eval_data[eval_data.error > 0]

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
        sun = self.star.sun

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

                assert np.all(
                    error > 0
                ), "Can only combine data points, not upper limits"

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
        self.det = det
        self.cov = cov

        del self.data
        del self.ions

        if show:
            print()
            self.print_db(ind=4, fields=True)
        self.show = show

    def run(
        self,
        stars,
        offsets=None,
        optimize=True,
        **kwargs,
    ):
        """Solve for the specified stars"""
        stars = np.array(stars)
        sol = np.ndarray(stars.shape, dtype=[("index", "int"), ("offset", "f8")])

        # TODO - rewrite loop to do vector operations if possible
        for i in range(stars.shape[0]):
            sol[i, :]["index"] = stars[i, :]
            if offsets is None:
                assert optimize is True
                sol[i, :]["offset"] = 1.0e-4 / stars.shape[1]
                if stars.shape[1] == 1:
                    local_search = False
                else:
                    local_search = True
            else:
                sol[i, :]["offset"] = offsets[i]
                if optimize is True:
                    local_search = False
                else:
                    local_search = None

        fitness = get_fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            sol,
            self.ejecta,
            cdf=self.cdf,
            dst=self.dst,
            local_search=local_search,
            limit_solver=self.limit_solver,
            limit_solution=self.limit_solution,
            **kwargs,
        )
        return sol, fitness

    def text_result(
        self,
        n=20,
        *,
        n0=0,
        format="unicode",
        best=True,
        wide=12,
        show_index=True,
        _return_dbx=False,
    ):
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

        if not show_index:
            base_title = base_title[:-1]
            base_units = base_units[:-1]

        empty_title = [""] * len(base_title)

        if best:
            stars = self.sorted_stars
            fitness = self.sorted_fitness
        else:
            stars = self.unsorted_stars
            fitness = self.unsorted_fitness

        n1 = min(n0 + n, len(stars))

        dbidx = np.array(self.db_idx[stars[n0:n1]["index"]])
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
            line = base_title.copy()
            for word in title:
                if word in _title_translate:
                    word = _title_translate[word]
                    if format == "html":
                        word = _title_translate_html.get(word, word)
                elif np.all([w.isalpha() for w in word.split()]):
                    word = capwords(word)
                line.append(word)
            text.append(line)
            line = base_units.copy()
            for word in units:
                if word in ("-", ""):
                    word = ""
                else:
                    word = _unit_translate.get(word, word)
                    if format == "html":
                        word = _unit_translate_html.get(word, word)
                    word = f"({word})"
                line.append(word)
            text.append(line)

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
                index, offset = stars[i, j]
                db_idx = self.db_idx[index]
                db = self.db[db_idx]
                if db_idx != db_idx0 and not wide:
                    _short_head(db, j == 0 and self.sol_size > 1)
                    db_idx0 = db_idx
                dbindex = index - self.db_off[db_idx]
                data = db.fielddata[dbindex]
                line = list()
                if j == 0:
                    line.append(f"{fitness[i]:3.2f}")
                else:
                    line.append("")
                line.append(f"{np.log10(offset):7.2f}")
                if self.db_n > 1:
                    line.append(f"{self.db_lab[db_idx]}")
                if show_index:
                    line.append(f"{dbindex:6d}")
                if wide:
                    for k in range(nfield):
                        if k in fieldmap[db_idx]:
                            l = np.where(fieldmap[db_idx] == k)[0][0]
                            d = data[l]
                            f = db.fieldformats[l]
                            if isinstance(d, bytes):
                                d = d.decode()
                            line.append(f"{d:{f}}")
                        else:
                            line.append("")
                else:
                    for d, f in zip(data, db.fieldformats):
                        if isinstance(d, bytes):
                            d = d.decode()
                        line.append(f"{d:{f}}")
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

    def text_db(self, dbx=None, filename=False, fields=False):
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
                num = False
        lines = list()
        line = list()
        if not num:
            line.append("NR")
        line.append(f"{'DB':>{db_len}}")
        if filename:
            line.append("File")
        else:
            line.append("Name")
        lines.append(line)
        dbn_len = 0
        for i in dbx:
            db = self.db[i]
            line = list()
            if not num:
                line.append(f"{i:>2d}")
            line.append(f"{self.db_lab[i]:>{db_len}}")
            if filename:
                line.append(Path(db.filename).name)
            else:
                line.append(db.name)
            dbn_len = max(dbn_len, len(line[-1]))
            lines.append(line)
        if fields:
            lines[0][-1] += " " * (dbn_len - len(lines[0][-1]))
            lines[0].append("Fields")
            for j_, i in enumerate(dbx):
                db = self.db[i]
                j = j_ + 1
                fieldinfo = list()
                for f, u in zip(db.fieldnames, db.fieldunits):
                    if len(u) > 0:
                        fieldinfo.append(f"{f} [{u}]")
                    else:
                        fieldinfo.append(f)
                lines[j][-1] += " " * (dbn_len - len(lines[j][-1]))
                lines[j].append(", ".join(fieldinfo))
        return lines

    def format_db(self, ind=0, pad="", **kwargs):
        if pad is None:
            pad = ""
        pad = " " * ind + pad
        lines = self.text_db(**kwargs)
        return "\n".join([pad + "  ".join(line) for line in lines])

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

    def plot(
        self,
        num=0,
        fig=None,
        ax=None,
        fontsize=None,
        annosize=None,
        figsize=(10, 6),
        dpi=102,
        show_copyright=True,
        dist=None,
        xlim=None,
        ylim=None,
        data_size=3,
        save=None,
        save_format=None,
        return_plot_data=False,
        yscale=2,
        ynorm="Fe",
        range_det=False,  # adjust range to include detection thresholds
        range_lim=True,  # adjust range to incoude detection limits
        range_zmin=3,  # adjust range to consider abundances with at least that Z
        pad_abu=0.1,
        pad_det=0.05,
        multi=0,
        xlabel=None,
        ylabel=None,
    ):
        r"""
        Create a plot of the solution.

        Note:
            The legend as well as the star name and copyright string can be
            moved (dragged).

        Args:
            num (int): Number of solution, from the top (default: ``0``).
            yscale (int): select the y-scale of the plot.  Numerical value
                identical to those used for the star data formats.
            ynorm (str): elements to use a norm for ``[X/Y]`` plots
                (``yscale=3``).
            multi (int): plot this many best solutions as grey lines
                (default: ``0``).  For ``multi=-1`` lines will be shaded
                according to relative data point probability based on
                :math:`{\chi^2}`; and
                assuming multi-dimensional Gaussian error.
            save (str): filename to save plot.
            range_det (bool): adjust range to include detection thresholds
                (default: ``False``)
            range_lim: (bool) adjust range to include detection limits
                (default: ``True``)
            range_zmin (int): minimum Z to consider for determining y range
                (default: ``3``)
            pad_abu (float): fraction of plot range to use at bottom/top
                boundary (default: ``0.1``).
            pad_det (float): fraction of plot range to pad detection thresholds
                (default: ``0.05``).
            figsize (tuple): dimensions of figure in inches
                (default: ``(10, 6)``).
            dpi (int): resolution of image (default: ``102``).
            xlim (tuple): overwrite x range (low, high).
            ylim (tuple): overwrite y range (low, high).
            data_size (int): Size of data lines and symbols (default: ``3``).
            fontsize (int): size used for axis labels (default: ``12``).
            annosize (str): size used for element symbols (default: ``small``).
            dist (float):  distance of labels from data points.
            fig (:class:``matplotlib.figure.Figure``): figure object to use as
                canvas, otherwise as new figure is created.
            ax (:class:``matplotlib.axis.Axis``): axis objects to use for
                drawing, otherwise axis and parent figure are created as needed.
            xlabel (str): overwrite label for x-axis.
            ylabel (str): overwrite label for y-axis.
        """

        num = min(max(num, 0), len(self.sorted_stars) - 1)
        multi = min(max(multi, -1), len(self.sorted_stars))

        bestsol = self.sorted_stars[num]
        indices = bestsol["index"]
        offsets = bestsol["offset"]

        if ax is None:
            if fig is None:
                fig, ax = plt.subplots(
                    figsize=figsize,
                    dpi=dpi,
                    facecolor="white",
                    edgecolor="white",
                )
            else:
                ax = fig.add_subplot(111)

        if annosize is None:
            if fontsize is not None:
                annosize = fontsize
            else:
                annosize = "small"

        if fontsize is None:
            fontsize = 12

        convert = Convert_from_5(self, yscale=yscale, ynorm=ynorm, ylabel=ylabel)

        if xlabel is None:
            xlabel = "Element charge number"
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(convert.ylabel, fontsize=fontsize)
        ax.set_yscale(convert.plot_scale)

        zlist_db = np.array([ion.Z for ion in self.list_db])
        zlist_comb = np.array([ion.Z for ion in self.list_comb])

        # First, sum up the values for all the stars in the solution
        summed = np.sum((self.full_abudata[:, indices] + 1.0e-99) * offsets, axis=1)

        # Transform the list of matched elements to an index for the db
        index_t = np.in1d(zlist_db, zlist_comb, assume_unique=True)

        # logsun_full = np.log10(self.sun_full)
        # logsun_star = np.log10(self.sun_star)

        # The star data points
        y_star = convert.star_abundance
        y_star_det = convert.star_abundance_det

        x_star = ufunc_Z(self.eval_data.element)

        y_star_uco = np.abs(self.eval_data.error)
        y_star_cov = np.sqrt(np.sum(self.eval_data.covariance**2, axis=1))
        y_star_err = np.sqrt(y_star_cov**2 + y_star_uco**2)

        y_star_err = np.tile(y_star_err, (2, 1))
        y_star_uco = np.tile(y_star_uco, (2, 1))

        # set asymmetric error for upper limits
        up_lims = self.uplim_index_star
        y_star_err[1, up_lims] = 0

        cov_sel = (y_star_cov > 0) & ~up_lims

        # convert log errors into coordinate_errors
        y_star_err = convert.err(y_star_err)
        y_star_uco = convert.err(y_star_uco)

        leg_starname(ax, self.star.name)

        if show_copyright:
            leg_copyright(ax)

        # The pattern found by the algorithm
        y_a = convert(summed[index_t], ref="full", scale="lin")
        x_a = np.array(zlist_comb)

        x_ga = np.array(zlist_comb)

        # Move the x position for combined things
        x_star = x_star.astype("f8")
        x_a = x_a.astype("f8")
        x_ga = x_ga.astype("f8")
        for row in self.combine:
            if len(row) > 0:
                for x in [x_star, x_a, x_ga]:

                    # Get the first in the combined elements
                    ind = np.where(x == row[0])[0]
                    if len(ind) > 0:

                        # Make the new position the average of the first
                        # consecutives
                        consec = np.array_split(row, np.where(np.diff(row) != 1)[0] + 1)
                        x[ind[0]] = np.mean(consec[0])

        # superimpose alternate solutions by weight
        if multi < 0:
            w = self.W()
            wm = np.max(w)
            for i in range(self.db_size):
                if i == num:
                    continue
                alpha = w[i] / wm
                if alpha * 1024 < 1:
                    continue
                sol = self.unsorted_stars[i]
                sol_indices = sol["index"]
                sol_offsets = sol["offset"]
                sol_summed = np.sum(
                    (self.full_abudata[:, sol_indices] + 1.0e-99) * sol_offsets, axis=1
                )
                y_ga = convert(sol_summed[index_t], ref="full", scale="lin")
                ax.plot(
                    x_ga,
                    y_ga,
                    linestyle="-",
                    lw=2,
                    zorder=-10,
                    color="#cfcfcf",
                    alpha=alpha,
                )

        # plot alternative solutions
        for i in range(multi):
            if i == num:
                continue
            sol = self.sorted_stars[i]
            sol_indices = sol["index"]
            sol_offsets = sol["offset"]
            sol_summed = np.sum(
                (self.full_abudata[:, sol_indices] + 1.0e-99) * sol_offsets, axis=1
            )
            y_ga = convert(sol_summed[index_t], ref="full", scale="lin")
            ax.plot(
                x_ga,
                y_ga,
                linestyle="-",
                lw=2,
                zorder=-10,
                color="#cfcfcf3f",
            )

        # Components of the solution
        lines = ["--", "-.", ":", "-"]
        linecycler = cycle(lines)

        labels = list()
        texlabels = list()
        lines_used = list()

        for i, (offset, index) in enumerate(zip(offsets, indices)):
            raw = list()
            parameters = list()
            db_idx = self.db_idx[index]
            db = self.db[db_idx]
            dbindex = index - self.db_off[db_idx]
            if self.db_n > 1:
                db_name = self.db_lab[db_idx]
                try:
                    int(db_name)
                    db_name = f"DB {db_name}"
                except:
                    pass
                parameters.append(db_name)

            for j in range(db.nfield):
                if db.fieldflags[j] != StarDB.Flags.parameter:
                    continue
                value = db.fielddata[dbindex][j]
                if isinstance(value, bytes):
                    value = value.decode()
                unit = db.fieldunits[j]
                name = db.fieldnames[j]
                form = db.fieldformats[j]
                if unit == "-":
                    unit = ""
                raw.append(f"{value:{form}} {unit}".strip())

                if name in _plot_title_formatters:
                    value = _plot_title_formatters[name](value, form, unit)
                elif unit in _plot_unit_formatters:
                    value = _plot_unit_formatters[unit](value, form)
                else:
                    value = f"{value:{form}}"
                    unit = _plot_unit_translate.get(unit, unit)
                    if unit not in (
                        "",
                        "-",
                    ):
                        value = f"{value} {unit}"
                parameters.append(value)

            texlabels.append(", ".join(parameters))
            if self.db_n > 1:
                key = f"{db_name}, {dbindex}"
            else:
                key = f"{dbindex}"
            labels.append(f"{key}: " + ", ".join(raw))

            # Plot components if there are more than one
            if len(indices) > 1:
                if multi == 0:
                    y_ga = convert(
                        self.full_abudata[index_t, index] * offset,
                        ref="full",
                        scale="lin",
                    )
                    lines_used.append(next(linecycler))
                    ax.plot(
                        x_ga,
                        y_ga,
                        linestyle=lines_used[-1],
                        lw=1.3,
                        color=colorsys.hsv_to_rgb(i * 360 / len(indices), 1, 1),
                        # label=texlabels[i],
                    )
                else:
                    lines_used.append("none")

        # Change the label of the summed line based on how many components there are
        chi = rf"($\chi^2={self.sorted_fitness[num]:0.2f})$"
        if len(indices) > 1:
            sumlabel = f"Sum {chi}"
        else:
            sumlabel = f"{texlabels[0]} {chi}"

        # Hidden lines for legend purposes
        ax.plot(
            [None],
            [None],
            marker="^",
            markerfacecolor="red",
            markeredgecolor="red",
            color="g",
            lw=2.0,
            label=sumlabel,
        )

        if len(indices) > 1:
            for i, ls in enumerate(lines_used):
                if multi == 0:
                    ax.plot(
                        [None],
                        [None],
                        linestyle=ls,
                        lw=1.3,
                        color=colorsys.hsv_to_rgb(i * 360 / len(indices), 1, 1),
                        label=texlabels[i],
                    )
                else:
                    ax.plot(
                        [None],
                        [None],
                        linestyle="none",
                        lw=0,
                        color="none",
                        label=texlabels[i],
                    )

        if multi != 0 and not (multi == 1 and num == 0):
            if multi > 0:
                label = f"{multi:d} alternative solutions"
            else:
                label = "weighed alternative solutions"
            ax.plot(
                [None],
                [None],
                linestyle="-",
                lw=2,
                zorder=-10,
                label=label,
                color="#cfcfcfcf",
            )

        # Green line
        ax.plot(
            x_a,
            y_a,
            color="g",
            lw=2.0,
        )

        # Red triangles
        ax.plot(
            x_a[~self.lolim_index_all],
            y_a[~self.lolim_index_all],
            marker="^",
            markerfacecolor="red",
            markeredgecolor="red",
            lw=0,
        )

        # Hollow triangles
        ax.plot(
            x_a[self.lolim_index_all],
            y_a[self.lolim_index_all],
            marker="^",
            markerfacecolor="white",
            markeredgecolor="red",
            lw=0,
        )

        if dist is None:
            # Calculate bounding box
            txt = fig.text(0, 0, "Mg", fontsize=annosize)
            renderer = fig.canvas.get_renderer()
            bbox = txt.get_window_extent(renderer)
            txt.remove()
            dist = bbox.height

        # Plot limits
        if ylim is None:
            ii = x_star >= range_zmin
            if not range_lim:
                ii &= y_star_err[1, :] > 0.0
            ylim = np.array(
                [
                    np.min(y_star[ii] - y_star_err[0, ii]),
                    np.max(y_star[ii] + y_star_err[1, ii]),
                ]
            )
            if convert.plot_scale == "linear":
                dylim = ylim[1] - ylim[0]
                ylim += dylim * np.array([-1, 1]) * pad_abu
            else:
                dylim = np.log10(ylim[1]) - np.log10(ylim[0])
                ylim *= 10 ** (np.array([-1, 1]) * pad_abu * dylim)

            if range_det:
                y_min_det = np.array([y for y in y_star_det[ii] if y > convert.det_lim])
                if len(y_min_det) > 0:
                    if convert.plot_scale == "linear":
                        ylim[0] = min(ylim[0], min(y_min_det) - pad_det * dylim)
                    else:
                        ylim[0] = min(
                            ylim[0], min(y_min_det) * 10 ** (-pad_det * dylim)
                        )

        if xlim is None:
            xlim = (self.z_min - 0.5, self.z_max + 0.5)

        # Calculate number of pixels per data
        dpi = fig.get_dpi()
        height_inch = figsize[1]
        height_pixel = height_inch * dpi
        if convert.plot_scale == "linear":
            data_scale = height_pixel / (ylim[1] - ylim[0])
            space = fontsize / data_scale
            gap_size = 3.5 * space
        else:
            data_scale = height_pixel / (np.log10(ylim[1]) - np.log10(ylim[0]))
            space = fontsize / data_scale
            gap_size = 10 ** (3.5 * space)

        anno = np.copy(y_a)

        # This is the ind corresponding to elements in the star data (the error bar points)
        for ind, (z, zor) in enumerate(zip(x_a, zlist_comb)):
            # We use x_a instead of zlist_comb because it has the modified Z
            # for the combined elements

            y = y_a[ind]

            # determine whether to plot above or below
            if ind == 0:
                above = y_a[ind + 1] < y
            elif ind == len(zlist_comb) - 1:
                above = y_a[ind - 1] < y
            else:
                above = 2 * y >= y_a[ind - 1] + y_a[ind + 1]

            if above:
                loc = (0, 0.7 * dist)
            else:
                loc = (0, -dist)

            if z in x_star:
                star_ind = np.where(x_star == z)[0][0]
                ys = y_star[star_ind]
                ye = y_star_err[:, star_ind]

                if above:
                    y_up = ys + ye[1]
                    if convert.plot_scale == "linear":
                        upper = y + gap_size
                    else:
                        upper = y * gap_size
                    if (y < y_up) and (upper > ys - ye[0]):
                        anno[ind] = y_up
                else:
                    y_lo = ys - ye[0]
                    if convert.plot_scale == "linear":
                        lower = y - gap_size
                    else:
                        lower = y / gap_size
                    if (y > y_lo) and (lower < ys + ye[1]):
                        anno[ind] = y_lo

            # Make a special label for combined things
            label = None
            for row, group in enumerate(self.combine):
                if len(group) > 0:
                    if zor == group[0]:
                        strings = [str(I(z)) for z in self.combine[row]]
                        label = "+".join(strings)
            if label is None:
                label = I(int(z)).element_symbol()

            annotation = ax.annotate(
                label,
                xy=(z, anno[ind]),
                xycoords="data",
                xytext=loc,
                textcoords="offset points",
                size=annosize,
                ha="center",
                va="center",
                clip_on=True,
                annotation_clip=False,
            )

            annotation.draggable(True)

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
            if y < -20:
                continue
            ax.plot(
                x + 0.4 * np.array([-1, 1]),
                np.array([y, y]),
                ls="-",
                lw=data_size,
                color="#0000FF3f",
            )

        # Show correlated errors
        ax.errorbar(
            np.array(x_star)[cov_sel],
            y_star[cov_sel],
            yerr=y_star_uco[:, cov_sel],
            ls="None",
            marker=None,
            ms=0,
            lw=2 * data_size,
            capsize=0,
            color="#0000ff3f",
            mfc=(0, 0, 0),
            uplims=False,
            zorder=-1,
        )

        # Plot for the excluded data points
        ax.errorbar(
            np.array(x_star)[self.exclude_index],
            y_star[self.exclude_index],
            yerr=y_star_err[:, self.exclude_index],
            ls="None",
            marker="o",
            ms=data_size * 1.5,
            color=(0, 0, 0),
            mfc=(1, 1, 1),
            capsize=data_size,
            uplims=up_lims[self.exclude_index],
        )

        # Plot for the data points
        ax.errorbar(
            np.array(x_star)[~self.exclude_index],
            y_star[~self.exclude_index],
            yerr=y_star_err[:, ~self.exclude_index],
            ls="None",
            marker="o",
            ms=0,
            capsize=data_size,
            color=(0, 0, 0),
            mfc=(0, 0, 0),
            uplims=up_lims[~self.exclude_index],
        )

        # Make some dots for points that aren't uplims
        ii = np.where(~self.exclude_index)[0]
        ij = np.where(up_lims[ii] == 0)[0]
        ii = ii[ij]
        ax.scatter(
            np.array(x_star)[ii],
            y_star[ii],
            marker="o",
            s=10 * data_size,
            color=(0, 0, 0),
        )

        ax.set_ybound(*ylim)
        ax.set_xbound(*xlim)

        # ax.ticklabel_format(style='plain')
        ax.tick_params(axis="both", which="major", labelsize=fontsize)

        ax.xaxis.set_major_formatter(IntFormatter())
        if convert.plot_scale == "linear":
            ax.yaxis.set_major_formatter(IntFormatter())

        fig.tight_layout()
        fig.show()

        if save is True:
            save = (
                self.star.name
                + "."
                + "-".join([str(index) for index in indices])
                + "."
                + str(num)
            )
            save = Path.cwd().parent / "plots" / (save + ".pdf")
        if isinstance(save, (str, Path, BytesIO)):
            fig.savefig(save, format=save_format)

        if return_plot_data:
            return labels, (x_a, y_a)

    def plot_error_matrix(self, num=0, zoom=False, nlab=9):
        stars = self.sorted_stars[[num]]
        sol, fitness = self.run(
            stars=stars["index"],
            offsets=stars["offset"],
            fixed_offsets=False,
            optimize=False,
            return_matrix=True,
        )
        m = fitness[0]

        mx = m.copy()
        ii = np.arange(mx.shape[0])
        mx[ii, ii] = 0
        matrix = np.any(mx != 0)

        if zoom is not False:
            y = np.sign(m) * (np.log10(np.abs(m) + 1 / zoom) + np.log10(zoom))
            scale = " (scaled around zero)"
        else:
            y = m
            scale = ""
        vmag = np.max(np.abs(y))

        data = self.eval_data[~self.exclude_index]
        ions = data.element
        ticks = np.arange(len(ions)) + 0.5
        labels = [i.Name() for i in ions]

        label = rf"error contribution to $\chi^2${scale}"

        fig, ax = plt.subplots()
        if matrix:
            cm = ax.pcolormesh(y, vmin=-vmag, vmax=vmag, cmap="bwr")
            cb = fig.colorbar(cm)
            cb.set_label(label)
            if zoom is not False:
                x = ((np.arange(nlab) / (nlab - 1)) * 2 - 1) * vmag
                cb.set_ticks(x)
                y = np.sign(x) * (10 ** (np.abs(x) - np.log10(zoom)) - 1 / zoom)
                cb.set_ticklabels([f"{i:5.3f}" for i in y])
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
        else:
            values = np.diag(y)
            l = ax.plot([None])
            c = l[0].get_color()
            for x, y in zip(ticks, values):
                mfc = None if y > 0 else "none"
                ax.plot(x, y, ls="none", color=c, marker="o", mfc=mfc)
            ax.axhline(0, color="#cfcfcf", ls="-", zorder=-10)
            ax.set_ylabel(label)

        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)

        fig.tight_layout()

    def plot_star_matrix(self, zoom=False, nlab=9, compress=True):

        m = get_complete_matrix(self.star, self.cdf)
        if zoom is not False:
            y = np.sign(m) * (np.log10(np.abs(m) + 1 / zoom) + np.log10(zoom))
            scale = " (scaled around zero)"
        else:
            y = m
            scale = ""

        ions = self.star.element_abundances.element
        if compress:
            ii = np.any(y != 0.0, axis=0)
            y = y[ii][:, ii]
            ions = ions[ii]

        vmag = np.max(np.abs(y))

        fig, ax = plt.subplots()
        cm = ax.pcolormesh(y, vmin=-vmag, vmax=vmag, cmap="bwr")
        cb = fig.colorbar(cm)
        cb.set_label(f"matrix elements{scale}")
        if zoom is not False:
            x = ((np.arange(nlab) / (nlab - 1)) * 2 - 1) * vmag
            cb.set_ticks(x)
            y = np.sign(x) * (10 ** (np.abs(x) - np.log10(zoom)) - 1 / zoom)
            cb.set_ticklabels([f"{i:5.3f}" for i in y])

        ticks = np.arange(len(ions)) + 0.5
        labels = [i.Name() for i in ions]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

        fig.tight_layout()

    def plot_star_inverse(self, zoom=False, nlab=9, compress=True):

        m = get_complete_inverse(self.star, self.cdf)
        if zoom is not False:
            y = np.sign(m) * (np.log10(np.abs(m) + 1 / zoom) + np.log10(zoom))
            scale = " (scaled around zero)"
        else:
            y = m
            scale = ""

        ions = self.star.element_abundances.element
        if compress:
            ii = np.any(y != 0.0, axis=0)
            y = y[ii][:, ii]
            ions = ions[ii]

        vmag = np.max(np.abs(y))

        fig, ax = plt.subplots()
        cm = ax.pcolormesh(y, vmin=-vmag, vmax=vmag, cmap="bwr")
        cb = fig.colorbar(cm)
        cb.set_label(f"inverse matrix elements{scale}")
        if zoom is not False:
            x = ((np.arange(nlab) / (nlab - 1)) * 2 - 1) * vmag
            cb.set_ticks(x)
            y = np.sign(x) * (10 ** (np.abs(x) - np.log10(zoom)) - 1 / zoom)
            cb.set_ticklabels([f"{i:5.3f}" for i in y])

        ticks = np.arange(len(ions)) + 0.5
        labels = [i.Name() for i in ions]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

        fig.tight_layout()

    def W(self):
        k = self.fit_size - 1
        k2 = k * 0.5
        chi2 = self.unsorted_fitness * k
        return chi2 ** (k2 - 1) * np.exp(-0.5 * chi2) * np.sqrt(0.5) ** k / gamma(k2)


class Convert_from_5(object):
    def __init__(self, starfit, yscale=2, ynorm=None, ylabel=None):
        self.yscale = yscale
        self.ynorm = ynorm

        self.star_log_abundance = starfit.eval_data.abundance
        self.star_log_abundance_det = starfit.eval_data.detection

        if self.yscale == 1:
            if starfit.star.data_format == 1:
                self.log_H = np.log10(starfit.star.norm_element)
            else:
                self.log_H = starfit.star.element_abundances.abundance[0]
        if self.yscale in (2, 3, 6):
            self.log_sun_full = np.log10(starfit.sun_full)
            self.log_sun_star = np.log10(starfit.sun_star)
        if self.yscale in (4, 7):
            self.ynorm = "Si"
        if self.yscale == 6:
            self.ynorm = "H"
        if starfit.star.data_format == 3:
            if self.ynorm is None:
                self.ynorm = starfit.star.norm_element.Name()
        if self.yscale in (3, 4, 6, 7):
            i = np.where(starfit.star.element_abundances.element == I(self.ynorm))[0]
            if len(i) == 0:
                raise AttributeError(f"Norm Element {self.ynorm} not available")
            i = i[0]
            self.log_norm_star = starfit.star.element_abundances.abundance[i]
            self.log_norm_sun = np.log10(starfit.star.sun.Y(self.ynorm))

        self.star_abundance = self.__call__(self.star_log_abundance, "star")
        self.star_abundance_det = self.__call__(self.star_log_abundance_det, "star")

        if self.yscale in (0, 7):
            self.plot_scale = "log"
            self.det_lim = 1e-20
        else:
            self.plot_scale = "linear"
            self.det_lim = -20

        self._ylabel = ylabel

    def __call__(self, y, ref=None, scale="log"):
        if scale == "lin":
            y = np.log10(y + 1.0e-99)
        if self.yscale == 5:
            return y
        if self.yscale == 0:
            return 10**y
        if self.yscale == 1:
            return y + self.log_H + 12
        if self.yscale in (2, 3, 6):
            if ref == "full":
                y = y - self.log_sun_full
            elif ref == "star":
                y = y - self.log_sun_star
        if self.yscale == 2:
            return y
        if self.yscale in (3, 6):
            return y - self.log_norm_star + self.log_norm_sun
        if self.yscale in (4, 7):
            y = y - self.log_norm_star + 6
            if self.yscale == 4:
                return y
            else:
                return 10**y

    _labels = {
        0: "Abundace",
        1: "Abundace (log epsilon)",
        2: "Logarithm of abundance relative to sun",
        3: "[X/{norm}]",
        4: "Logarithm of abundance ({norm} = 10$^{{6}}$ atoms)",
        5: "Logarithm of abundance (mol$\\,$/$\\,$g)",
        6: "[X/{norm}]",
        7: "Abundance ({norm} = 10$^{{6}}$ atoms)",
    }

    @property
    def ylabel(self):
        if self._ylabel is not None:
            return self._ylabel
        return self._labels[self.yscale].format(norm=self.ynorm)

    def err(self, err):
        err = err + self.star_log_abundance
        err = self.__call__(err, "star")
        err = err - self.star_abundance
        return err
