"""Genetic algorithm"""

import time

import numpy as np

from .autils.human import time2human
from .autils.logged import Logged
from .starfit import Results, _fitness
from .starplot import fitplot


class Ga(Results, Logged):
    """Genetic Algorithm."""

    def __init__(
        self,
        filename,
        db="znuc.S4.star.el.y.stardb.gz",
        z_min=1,
        z_max=30,
        z_exclude=[],
        z_lolim=[],
        combine=[[]],
        gen=10,
        time_limit=20,
        pop_size=200,
        sol_size=None,
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
        max_pop=2**13,
        cover=None,
    ):

        t_total = 0
        self.sol_size = sol_size
        self.silent = silent
        super().__init__()

        self.rng = np.random.default_rng(seed)

        self._setup(
            filename=filename,
            database=db,
            combine=combine,
            z_exclude=z_exclude,
            z_min=z_min,
            z_max=z_max,
            upper_lim=upper_lim,
            z_lolim=z_lolim,
        )

        if sol_size is None:
            if cover is None:
                cover = True
            if self.db_n > 1 and cover:
                sol_size = self.db_n
            else:
                sol_size = 2
        if cover is None:
            cover = False

        # If offsets is fixed, there is no offset mutation
        if fixed_offsets:
            mut_rate_offset = 0
            local_search = False

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
        self.cover = cover

        self.logger.info(f"Time limit: {time2human(time_limit)}")

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
        self.s = self._populate()

        # Generate initial fitness values
        self.f = _fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            self.s,
            fixed_offsets=fixed_offsets,
            ejecta=self.ejecta,
            cdf=cdf,
            ls=self.local_search,
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
        shuffle = self.rng.permutation(ft.shape[0])
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
        if self.cover:
            mut_db = self.rng.integers(self.db_n, size=w.shape)
            mutants = self.db_off[mut_db] + self.rng.integers(self.db_num[mut_db])
        else:
            mutants = self.rng.integers(self.db_size, size=w.shape)

        w["index"] = np.choose(
            self.rng.uniform(size=w.shape) < self.mut_rate_index,
            [
                w["index"],
                mutants,
            ],
        )
        if not self.local_search:
            w["offset"] = np.choose(
                self.rng.uniform(size=w.shape) < self.mut_rate_offset,
                [
                    w["offset"],
                    np.exp(self.rng.normal(self.mut_offset_magnitude, size=w.shape))
                    * w["offset"],
                ],
            )

        # Single point crossover
        w = w.reshape(2, -1, self.sol_size)
        crossind = self.rng.integers(self.sol_size, size=w.shape[1])

        # Make children and inverse-children
        c = np.ndarray(
            (w.shape[0] * w.shape[1], self.sol_size),
            dtype=[("index", np.int64), ("offset", np.float64)],
        )
        for i in range(0, w.shape[1]):
            cut = crossind[i]
            c[2 * i, :cut] = w[0, i, :cut]
            c[2 * i, cut:] = w[1, i, cut:]
            c[2 * i + 1, :cut] = w[1, i, :cut]
            c[2 * i + 1, cut:] = w[0, i, cut:]

        f = np.copy(f)
        fc = _fitness(
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

        # Sort entries by index
        o = np.take_along_axis(o, np.argsort(o["index"], -1), -1)

        # Eliminate genes with replications
        idx = o["index"]
        mask = np.all(idx[:, :-1] != idx[:, 1:], axis=-1)
        o = o[mask]
        f = f[mask]

        # If requested, elminate solutions that do not span all DBs.
        if self.cover:
            idb = self.db_idx[o["index"]]
            mask = np.all(idb[:, 1:] <= idb[:, :-1] + 1, axis=-1)
            mask &= idb[:, 0] == 0
            mask &= idb[:, -1] == self.db_n - 1
            o = o[mask]
            f = f[mask]
            # n = np.count_nonzero(~mask)
            # if n > 0:
            #     self.logger.info(f"cover: eliminating {n} candidates, {f.shape[0]} remaining.")
            assert f.shape[0] > 0, "cover: no data left candidate elimination"

        # Duplicate elimination
        if self.local_search or self.fixed_offsets:
            idx = o["index"]
            ind = np.lexsort(idx.T[::-1])
            sel = ind[
                np.concatenate(
                    (
                        [True],
                        np.any(
                            idx[ind[1:]] != idx[ind[:-1]],
                            axis=1,
                        ),
                    )
                )
            ]
        else:
            sel = np.unique(f, return_index=True)[1]
        # n = f.shape[0] - len(sel)
        o = o[sel]
        f = f[sel]
        # self.logger.info(f"eliminating {n} duplicates, {f.shape[0]} remaining.")
        assert f.shape[0] > 0, "no data left after elimiation of duplicates"

        # Selection
        # Order by fitness
        sort = np.argsort(f)
        o = o[sort]
        f = f[sort]

        # Discard some, but not the elite ones
        elite_point = int(self.pop_size * self.frac_elite)
        if o.shape[0] > self.pop_size:
            discard_index = self.rng.choice(
                np.arange(elite_point, o.shape[0]),
                size=(o.shape[0] - self.pop_size,),
                replace=False,
            )
            s = np.delete(o, discard_index, axis=0)
            f = np.delete(f, discard_index, axis=0)
        else:
            self.max_pop = self.pop_size
            s = o

        # Trim solution vector to even size
        if len(f) % 2 != 0:
            f = f[:-1]
            s = s[:-1]

        self.s = s
        self.f = f

    def _populate(self):
        """Create an initial population"""

        assert (
            self.db_size >= self.sol_size
        ), f"too few db entries ({self.db_size}) for solution size ({self.sol_size})."

        s = np.ndarray(
            (self.pop_size, self.sol_size),
            dtype=[("index", np.int64), ("offset", np.float64)],
        )

        if self.cover:
            assert np.all(
                self.db_num > (self.sol_size - self.db_n + 1)
            ), "at least one db has too few elemets for requested solution size"
            idb = self.rng.permuted(
                np.tile(
                    np.arange(self.db_n, dtype=np.int64),
                    (self.pop_size, 1),
                ),
                axis=-1,
            )[:, : self.sol_size]
            if self.db_n < self.sol_size:
                idb = np.concatenate(
                    (
                        idb,
                        self.rng.integers(
                            self.db_n,
                            size=(
                                self.pop_size,
                                self.sol_size - self.db_n,
                            ),
                        ),
                    ),
                    axis=-1,
                )
            idx = self.db_off[idb] + self.rng.integers(self.db_num[idb])

            # replace duplicates
            while True:
                idx_sorted = np.sort(idx, axis=-1)
                ii = np.any(idx_sorted[:, 1:] == idx_sorted[:, :-1], axis=-1)
                if not np.any(ii):
                    break
                idx[ii] = self.db_off[idb[ii]] + self.rng.integers(self.db_num[idb[ii]])

            s["index"] = idx
        else:
            # ensure no duplicates
            n = 0
            while n < self.pop_size:
                m = min(self.db_size // self.sol_size, self.pop_size - n)
                s["index"][n : n + m, :] = self.rng.permutation(self.db_size)[
                    : m * self.sol_size
                ].reshape(-1, self.sol_size)
                n += m
        s["offset"] = (
            self.rng.uniform(size=(self.pop_size, self.sol_size))
            * (0.9 / self.sol_size)
            + 1.0e-14
        )

        return s

    def plot_fitness(self):
        # Fitness over time plot
        fitplot(
            starname=self.star.name,
            generations=self.gen,
            popsize=self.pop_size,
            genesize=self.sol_size,
            times=self.times,
            history=self.history,
        )
