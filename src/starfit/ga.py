"""Genetic algorithm"""

import time

import numpy as np

from .autils.logged import Logged
from .autils.time2human import time2human
from .starfit import Results


class Ga(Results, Logged):
    """Genetic Algorithm."""

    def __init__(
        self,
        filename,
        db="znuc.S4.star.el.y.stardb.gz",
        z_max=30,
        z_exclude=[],
        z_lolim=[],
        combine=[[]],
        gen=10,
        time_limit=20,
        pop_size=200,
        sol_size=2,
        tour_size=2,
        frac_mating_pool=1,
        frac_elite=0.5,
        mut_rate_index=0.2,
        mut_rate_offset=0.1,
        fixed_offsets=False,
        mut_offset_magnitude=1,
        local_search=True,
        upper_lim=True,
        cdf=True,
        seed=None,
        silent=False,
        interpolate=False,
        max_pop=10000,
    ):

        t_total = 0
        self.silent = silent
        Results.__init__(self, "GA")

        if seed is not None:
            np.random.seed(seed)

        self._setup(
            filename,
            db,
            combine,
            z_exclude,
            z_max,
            upper_lim,
            z_lolim,
            interpolate,
        )

        self.sol_size = sol_size
        self.pop_size = pop_size
        self.frac_elite = frac_elite
        self.local_search = local_search
        self.fixed_offsets = fixed_offsets
        self.time_limit = time_limit
        self.tour_size = tour_size
        self.mut_rate_index = mut_rate_index
        self.mut_rate_offset = mut_rate_offset
        self.mut_offset_magnitude = mut_offset_magnitude
        self.cdf = cdf
        self.max_pop = max_pop

        self.logger.info(f"Time limit: {time2human(time_limit)}")

        # If offsets is fixed, there is no offset mutation
        if fixed_offsets:
            mut_rate_offset = 0
            local_search = False

        # Turn all of the fractions into numbers
        # The tournament uses these numbers to determine the size of the output
        # If frac_mating_pool is unspecified, skip this part
        if frac_mating_pool is not None:
            if frac_mating_pool >= 0:
                # Remove the appropriate fraction of entries, rounding down so
                # that each tournament group has enough participants
                self.n_winners = int(frac_mating_pool * pop_size)
                self.n_winners -= self.n_winners % (tour_size)

            # frac_mating_pool is a fraction of the number of solutions
            else:
                raise Exception(
                    f"frac_mating_pool = {frac_mating_pool:f}, must be greater than 0"
                )
        else:
            self.n_winners = pop_size

        # Making sure that there are an even number of children for binary crossover
        if self.n_winners % 2 != 0:
            self.n_winners -= tour_size

        # Generate initial population
        self.s = self._populate(
            self.db.data.transpose().shape[1],
            pop_size,
            sol_size,
        )

        # Generate initial fitness values
        self.f = self._fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            self.s,
            fixed_offsets=fixed_offsets,
            ejecta=self.ejecta,
            cdf=cdf,
        )

        # Initial history points
        self.history["best"] += [np.min(self.f)]
        self.history["average"] += [np.mean(self.f)]
        self.history["worst"] += [np.max(self.f)]
        self.initsol = self.s[np.argmin(self.f)]
        self.times += [0]

        # Iterate
        i = 0

        while i <= gen or time_limit is not None:
            t_start = time.time()
            self._step()
            self.history["best"] += [self.f[0]]
            self.history["average"] += [np.mean(self.f)]
            self.history["worst"] += [self.f[-1]]

            t_total += time.time() - t_start
            self.times += [t_total]
            i += 1
            if time_limit is not None:
                gen = i
                if t_total >= time_limit:
                    break

        self.sorted_stars = self.s
        self.sorted_fitness = self.f

        self.gen = gen
        self.logger.info("Best fitness:")
        self.logger.info(f" {self.sorted_fitness[0]}")
        self.close_logger(timing="Finished in")

    def _step(self):
        """Perform one step of the genetic algorithm"""
        f = self.f
        s = self.s

        # Tournament
        t = np.repeat(s, self.tour_size, axis=0)
        ft = np.repeat(f, self.tour_size, axis=0)
        # Shuffle by generating a random permutation
        shuffle = np.random.permutation(ft.shape[0])
        t = t[shuffle]
        ft = ft[shuffle]
        # Turn the tournament pool into pairs (or whatever the tournament size is)
        # Keep rows according to n_winners only
        t = np.swapaxes(
            t.reshape(self.tour_size, -1, self.sol_size)[:, 0 : self.n_winners],
            1,
            2,
        )
        ft = ft.reshape(-1, self.tour_size)[0 : self.n_winners]
        # Find winners
        w = np.choose(np.argmin(ft, axis=1), t).transpose()
        # Mutate
        # Mutation is performed using a crossover-like method,
        # between the original array, and a fully mutant array
        w["index"] = np.choose(
            np.random.ranf(w.shape) < self.mut_rate_index,
            [
                w["index"],
                np.random.randint(1, self.db.data.transpose().shape[1], w.shape),
            ],
        )
        if not self.local_search:
            w["offset"] = np.choose(
                np.random.ranf(w.shape) < self.mut_rate_offset,
                [
                    w["offset"],
                    np.exp(np.random.normal(0, self.mut_offset_magnitude, w.shape))
                    * w["offset"],
                ],
            )
        # Singlepoint crossover
        w = w.reshape(2, -1, self.sol_size)
        crossind = np.random.randint(self.sol_size, size=w.shape[1])

        # Make children and inverse-children
        c = np.ndarray(
            (w.shape[0] * w.shape[1], self.sol_size),
            dtype=[("index", "int"), ("offset", "f8")],
        )
        for i in range(0, w.shape[1]):
            cut = crossind[i]
            c[2 * i, :cut] = w[0, i, :cut]
            c[2 * i, cut:] = w[1, i, cut:]
            c[2 * i + 1, :cut] = w[1, i, :cut]
            c[2 * i + 1, cut:] = w[0, i, cut:]

        f = np.copy(f)
        fc = self._fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            c,
            fixed_offsets=self.fixed_offsets,
            ejecta=self.ejecta,
            cdf=self.cdf,
            ls=self.local_search,
        )

        # Recombine
        o = np.concatenate((s, c), axis=0)
        f = np.concatenate((f, fc))

        # Duplicate elimination
        if self.local_search or self.fixed_offsets:
            oindex = o["index"]
            ind = np.lexsort(oindex.transpose())
            unq = ind[
                np.concatenate(
                    (
                        [True],
                        np.any(
                            oindex[ind[1:]] != oindex[ind[:-1]],
                            axis=1,
                        ),
                    )
                )
            ]
        else:
            unq = np.unique(f, return_index=True)[1]
        o = o[unq]
        f = f[unq]

        # Selection
        # Order by fitness
        sort = np.argsort(f)
        o = o[sort]
        f = f[sort]

        # Discard some, but not the elite ones
        elite_point = int(self.pop_size * self.frac_elite)
        if o.shape[0] > self.pop_size:
            discard_index = np.random.choice(
                np.arange(elite_point, o.shape[0]),
                o.shape[0] - self.pop_size,
                replace=False,
            )
            s = np.delete(o, discard_index, axis=0)
            f = np.delete(f, discard_index, axis=0)
        else:
            self.max_pop = self.pop_size
            s = o

        self.s = s
        self.f = f

        iev = self.likelihood(f)
        self.evidence = np.sum(iev)

        growth_condition = (
            iev[int(self.pop_size * self.frac_elite)] * self.n_comb()
            > (self.evidence / self.pop_size) ** 2
        )
        if growth_condition and (self.pop_size < self.max_pop):
            self._enlarge()

    def _enlarge(self):
        """Double the population size"""
        db_size = self.db.data.transpose().shape[1]
        old_pop_size = self.pop_size
        self.pop_size *= 2
        s = np.ndarray(
            (self.pop_size, self.sol_size), dtype=[("index", "int"), ("offset", "f8")]
        )
        f = np.ndarray(
            (self.pop_size),
        )
        s[:old_pop_size] = self.s
        s[old_pop_size:] = self._populate(db_size, old_pop_size, self.sol_size)
        f[:old_pop_size] = self.f
        f[old_pop_size:] = self._fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            s[old_pop_size:],
            fixed_offsets=self.fixed_offsets,
            ejecta=self.ejecta,
            cdf=self.cdf,
            ls=self.local_search,
        )

        # Order by fitness
        sort = np.argsort(f)
        s = s[sort]
        f = f[sort]

        self.s = s
        self.f = f

    @staticmethod
    def _populate(dbsize, pop_size, sol_size):
        """Create an initial population"""

        s = np.ndarray((pop_size, sol_size), dtype=[("index", "int"), ("offset", "f8")])
        s["index"] = np.random.randint(
            1,
            dbsize,
            (pop_size, sol_size),
        )
        s["offset"] = np.random.rand(pop_size, sol_size) * 0.1 + 1.0e-14

        return s
