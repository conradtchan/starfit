"""Results objects from the various algorithms"""

import math
import os

import numpy as np

from . import DATA_DIR, starplot
from .autils import abusets
from .autils.isotope import Ion
from .autils.logged import Logged
from .dbtrim import TrimDB as StarDB
from .fitness import solver
from .read import Star


class Results(Logged):
    """
    Object for running the various algorithms and a container for the results
    """

    def __init__(self, sol_type):
        self.sol_type = sol_type
        self.history = {"best": [], "average": [], "worst": []}
        self.initsol = None
        self.bestsol = None
        self.times = []
        self.gen = None
        self.pop_size = None
        if sol_type == "Single":
            self.sol_size = 1
        elif sol_type == "Double":
            self.sol_size = 2
        else:
            self.sol_size = None

        try:
            self.silent
        except:
            self.silent = False

        self.setup_logger(silent=self.silent)

    def _setup(
        self,
        filename,
        db,
        combine,
        z_exclude,
        z_max,
        upper_lim,
        z_lolim,
        interpolate=False,
    ):
        """Prepare the data for the solvers. Trims the databases and excludes
        elements.
        """
        if not os.path.isfile(db):
            _db = os.path.join(DATA_DIR, "db", db)
            if os.path.isfile(_db):
                db = _db
            else:
                raise IOError(f"file {db} not found")
        # Read a star
        star = Star(filename, silent=self.silent)
        self.star = star
        # Read the database
        # Note on reading db:
        # db.data[element,entry]
        # For element 1-83
        # Entry 1-16800
        db = StarDB(db, silent=self.silent)
        if interpolate:
            self.logger.info("Using interpolated database")
            db.inter_db(new_energy=interpolate[0], new_mixing=interpolate[1])
        self.db = db

        self.ejecta = db.fielddata["mass"] - db.fielddata["remnant"]

        self.combine = combine
        self.z_max = z_max
        self.upper_lim = upper_lim
        self.z_lolim = z_lolim
        self.z_exclude = z_exclude

        # Remove elements with Z > zmax
        mask_zmax = np.array(
            [ion.Z <= z_max for ion in star.element_abundances.element]
        )
        eval_data = star.element_abundances[mask_zmax]
        # eval_data = eval_data #The other one gets chopped and changed

        # Remove upper limit elements if upper limits is not enabled
        mask_uplim = np.array([error > 0 for error in eval_data.error])
        if not upper_lim:
            eval_data = eval_data[mask_uplim]

        # List of DB ions which will get chopped
        dbions = self.db.ions

        # The full set of abundance data from the database
        full_abudata = np.copy(db.data.transpose())
        full_ions = np.copy(db.ions)

        # List of every every element in the DB
        list_db = np.array([ion for ion in db.ions])

        # List of every element in the DB less than zmax
        list_zmax = list_db[np.where(np.array([ion.Z for ion in list_db]) <= z_max)]

        # Prepare the sun
        sun = abusets.SolAbu(
            name=os.path.join(DATA_DIR, "ref/sollo22.dat"),
            silent=self.silent,
        )
        # Transforms the uncombined element numbers for use in the sun
        index_sun = np.in1d(list_db, list_zmax, assume_unique=True)
        sunions = dbions[index_sun]
        sun_full = sun.Y(db.ions[index_sun])
        sun_star = sun.Y(eval_data.element[:])

        # Combine elements
        list_comb = list_zmax
        comb_abudata = full_abudata
        comb_ions = full_ions
        # Structure of list:
        #  [[6,7,8], [1,2]] will combine C+N+O and H+He
        if len(combine[0]) > 0:
            self.logger.info("Combining elements:")
            for group in combine:
                elements = [Ion(z) for z in group]
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

        self.logger.info(f"Matching {len(eval_data)} data points up to Z={z_max}:")

        self.logger.info("    " + " ".join([str(ion) for ion in eval_data.element]))

        self.logger.info(
            f"with {np.sum(eval_data.error < 0)} upper limits in the data:"
        )
        self.logger.info(
            "    "
            + " ".join([str(ion) for ion in eval_data.element[eval_data.error < 0]])
        )

        self.logger.info(f"and {np.sum(lolim_index_star)} lower limits in the models:")

        self.logger.info(
            "    " + " ".join([str(ion) for ion in eval_data.element[lolim_index_star]])
        )

        self.eval_data = eval_data
        self.trimmed_db = trimmed_db
        self.full_abudata = full_abudata
        self.list_db = list_db
        self.list_zmax = list_zmax
        self.list_comb = list_comb
        self.exclude_index = exclude_index
        self.lolim_index_all = lolim_index_all
        self.lolim_index_star = lolim_index_star
        self.sun = sun
        self.sun_full = sun_full
        self.sun_star = sun_star

    @staticmethod
    def _fitness(
        trimmed_db,
        eval_data,
        z_exclude_index,
        sol,
        ejecta=[],
        fixed_offsets=False,
        write=True,
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
        fitness /= eval_data.error.shape[0] - np.sum(z_exclude_index) - 2
        return fitness

    def n_comb(
        self,
    ):
        try:
            self._n_comb
        except:
            db_size = self.trimmed_db.shape[1]
            self._n_comb = np.product(
                np.arange(db_size - self.sol_size + 1, db_size + 1).astype("f8")
            ) / math.factorial(self.sol_size)
        return self._n_comb

    def ev_const(
        self,
    ):
        try:
            self._ev_const
        except:
            n_el = np.sum(np.invert(self.exclude_index))
            self._ev_const = (
                1
                / np.sqrt(2 * np.pi) ** (n_el - 2)
                * np.abs(np.product(self.eval_data["error"]))
            )
        return self._ev_const

    def likelihood(self, fitness):
        n_el = np.sum(np.invert(self.exclude_index))
        return self.ev_const() * np.exp(-0.5 * (n_el - 2) * fitness)

    def log_likelihood_unscaled(self, fitness):
        n_el = np.sum(np.invert(self.exclude_index))
        return -0.5 * (n_el - 2) * fitness

    def dbmf(
        self,
    ):
        '''Calculate the "database mass function"'''
        masses = self.db.grid()[2][0, 0, :]["mass"]
        n = len(masses)
        midpoints = np.ndarray((n + 1))
        midpoints[1:-1] = 0.5 * (masses[:-1] + masses[1:])
        midpoints[0] = masses[0] - (midpoints[1] - masses[0])
        midpoints[-1] = masses[-1] + (masses[-1] - midpoints[-2])
        width = midpoints[1:] - midpoints[:-1]

        density = 1 / width / n

        return masses, density

    def get_prior_weight(self, prior_name=None):
        """Get weight for prior (same as rejection rate)"""
        weight = np.ones(self.db.data.shape[0])

        if prior_name is not None:
            mass, dens = self.dbmf()
            uniform_rate = 1 / dens

            if prior_name == "uniform":
                distribution = np.ones_like(mass)
            elif prior_name == "log_uniform":
                distribution = 1 / mass
            elif prior_name == "salpeter":
                distribution = mass ** (-2.35)

            # Reshaping to map onto the database
            weight = weight.reshape(self.db.grid()[1].shape)
            weight[:, :, :] = uniform_rate[np.newaxis, np.newaxis, :] * distribution
            weight = weight.reshape(-1)

        # Normalize to sum to 1
        weight = weight / np.sum(weight)

        return weight

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

        fitness = self._fitness(
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

    def plot(self, **kwargs):
        """Call plotting routines to plot the best fit."""
        db = self.db

        bestsol = self.sorted_stars[0]
        self.labels, self.plotdata = starplot.abuplot(
            indices=bestsol["index"].tolist(),
            offsets=bestsol["offset"].tolist(),
            star=self.star,
            db=db,
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
        starplot.fitplot(
            sol_type=self.sol_type,
            starname=self.star.name,
            generations=self.gen,
            popsize=self.pop_size,
            genesize=self.sol_size,
            times=self.times,
            history=self.history,
        )

    def text_result(self, n=20, print_text=True):
        """Print data of best fit."""
        text = []
        text += [["Chi**2", "Index", "Mass", "Energy", "Mixing", "Dilution", "Ejecta"]]
        text += [["", "", "M_sun", "B", "log(M_He)", "", "log(M_sun)"]]
        db_field_index = [
            np.where(self.db.fieldnames == x)[0][0]
            for x in ("mass", "energy", "mixing", "remnant")
        ]
        for i in range(n):
            for j in range(self.sol_size):
                index, offset = self.sorted_stars[i][j]
                text += [
                    [
                        f"{self.sorted_fitness[i]:3.2f}",
                        f"{index:6d}",
                        f"{self.db.fielddata[index][db_field_index[0]]:3.2f}",
                        f"{self.db.fielddata[index][db_field_index[1]]:3.2f}",
                        "{:3.2f}".format(
                            math.log10(
                                self.db.fielddata[index][db_field_index[2]] + 1.0e-99
                            )
                        ),
                        f"1:{1 / offset:1.0f}",
                        "{:3.2f}".format(
                            self.db.fielddata[index][db_field_index[0]]
                            - self.db.fielddata[index][db_field_index[3]]
                        ),
                    ]
                ]
            if self.sol_size > 1:
                text += [""]

        return text

        string = "\n".join(text)
        if print_text:
            print(string)
        else:
            return string

    @staticmethod
    def textpad(s, n):
        nspace = n - len(s)
        return s + " " * nspace

    def __str__(self):
        text = self.text_result(10)
        string = ""
        for line in text:
            for word in line:
                string += self.textpad(word, 11)
            string += "\n"
        return string

    def __repr__(self):
        return f"{self.__class__.__name__}({self.star.name})"
