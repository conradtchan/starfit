"""Genetic algorithm"""

import sys
import time
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb

from .autils.human import time2human
from .autils.utils import is_iterable
from .fit import get_fitness
from .solgen._solgen import gen_slice
from .starfit import StarFit
from .starplot import leg_copyright, leg_info, leg_starname
from .utils import getch


class Ga(StarFit):
    """Genetic Algorithm."""

    def __init__(
        self,
        *args,
        gen=1000,
        time_limit=20,
        pop_size=200,
        sol_size=None,
        tour_size=2,
        frac_mating_pool=1,
        frac_elite=0.5,
        mut_rate_index=0.2,
        mut_rate_offset=0.1,
        mut_offset_magnitude=1,
        local_search=True,
        fixed_offsets=False,
        seed=None,
        max_pop=2**13,
        spread=None,
        group=None,
        pin=None,
        n_top=20,
        interactive=True,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        if self.show:
            return

        t_total = 0

        self.rng = np.random.default_rng(seed)

        if group is None:
            group = True
        if sol_size == 1:
            group = False

        self._setup_group(group)

        if sol_size is None:
            if spread is None:
                spread = True
            if self.group_n > 1 and spread:
                sol_size = self.group_n
            else:
                sol_size = 2
        if spread is None:
            spread = False

        # pinning
        if pin is None:
            pin = list()
        if isinstance(pin, (int, str)):
            pin = [pin]
        assert is_iterable(pin), "require list for {pin=}"
        for i, p in enumerate(pin):
            # convert db name string to group index if group of length 1
            if isinstance(p, str):
                for j, db in enumerate(self.db):
                    if p == db.name:
                        for k, g in enumerate(self.group):
                            if len(g) == 1 and g[0] == j:
                                pin[i] = p = k
                                break
                        if isinstance(p, int):
                            break
                    if isinstance(p, int):
                        break
            assert isinstance(p, int), f"require integer database index for pin ({p})"
            assert 0 <= p < self.group_n, f"pin needs to be valid group index ({p})"
        # assert len(pin) == len(set(pin)), f"require unique pin indices: {pin=}"
        assert len(pin) <= sol_size, f"can't pin more than {sol_size=} groups"
        pin_free = [p for p in range(self.group_n) if p not in pin]
        counter = Counter(pin)
        self.pin = np.array(list(counter.keys()), dtype=np.int64)
        self.pin_count = np.array(list(counter.values()), dtype=np.int64)
        self.pin_n = int(np.sum(self.pin_count))
        self.pin_free = np.array(pin_free, dtype=np.int64)
        self.pin_free_n = len(pin_free)
        self.pin_free_min = np.min(self.pin_free)
        self.pin_free_max = np.max(self.pin_free)

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
        self.max_pop = max_pop
        self.spread = spread
        self.n_top = n_top

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
        self.f = get_fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            self.s,
            fixed_offsets=fixed_offsets,
            ejecta=self.ejecta,
            cdf=self.cdf,
            limit_solution=self.limit_solution,
            limit_solver=self.limit_solver,
            local_search=self.local_search,
        )

        # Initial history points
        self.history = dict()
        self.history["best"] = [np.min(self.f)]
        self.history["mean"] = [np.mean(self.f)]
        self.history["median"] = [np.median(self.f)]
        self.history["worst"] = [np.max(self.f)]
        self.times = [0]

        if gen is None or gen <= 0:
            gen = 2**30

        # Iterate
        i = 0

        if not self.silent:
            self.print_header()
            self.print_empty_table()

        time_start = time.perf_counter()
        elapsed = 0
        self.elapsed = elapsed
        t_start = time.time()

        finish = False
        while True:
            self._step()
            self.history["best"] += [self.f[0]]
            self.history["mean"] += [np.mean(self.f)]
            self.history["median"] += [np.median(self.f)]
            self.history["worst"] += [self.f[-1]]

            t_new = time.time()
            t_total += t_new - t_start
            self.times += [t_total]
            t_start = t_new
            i += 1

            elapsed = time.perf_counter() - time_start
            if time_limit is not None:
                if elapsed >= time_limit:
                    finish = True
            if not self.silent and interactive:
                if len(getch()) > 0:
                    sys.stdout.write("\x1b[A")
                    finish = True
            if i == gen:
                finish = True
            if elapsed - self.elapsed > 1 or i == 1 or finish:
                self.n_solved = i
                if time_limit is None:
                    time_frac_done = 0
                else:
                    time_frac_done = elapsed / time_limit
                if gen > 2**29:
                    gen_frac_done = 0
                else:
                    gen_frac_done = i / gen
                self.frac_done = max(time_frac_done, gen_frac_done)

                try:
                    self.cycle_time = (
                        0.01 * (elapsed - self.elapsed) + 0.99 * self.cycle_time
                    )
                except:
                    self.cycle_time = 1.0e-10

                self.elapsed = elapsed

                self.gen = i

                if not self.silent:
                    self.print_update()

                if finish:
                    break

        self.sorted_stars = self.s
        self.sorted_fitness = self.f

        self.logger.info("Best fitness:")
        self.logger.info(f" {self.sorted_fitness[0]}")
        self.close_logger(timing="Finished in")

    def print_header(self):
        print("\n")
        print("=================================================================")
        print(
            f"Evolving Population of {self.pop_size:,d} combinations of {self.sol_size} stars"
        )
        if self.db_n == 1:
            print(f"Database: {self.db[0].name}")
        else:
            print("Databases:")
            self.print_db(ind=4)
        print(f"Star: {self.star.name}")
        print(
            f"Using {self.fit_size:d} elements ({np.sum(self.eval_data.error < 0):d} limits)"
        )
        if self.fixed_offsets:
            print("Offsets are ejecta mass")
        else:
            print("Offsets solved using Newton-Raphson")
        print("=================================================================")

        fmt_head_star = " Index  Offset"
        fmt_data_star = "{:>6d}  {:6.2f}"
        fmt_pad = " " * 5
        if self.db_n > 1:
            db_len = 2
            for l in self.db_lab:
                db_len = max(db_len, len(l))
            fmt_head_star = f"{'DB':>{db_len}} " + fmt_head_star
            fmt_data_star = f"{{:>{db_len}}} " + fmt_data_star
        fmt_head_star = fmt_pad + fmt_head_star
        fmt_data_star = fmt_pad + fmt_data_star
        self.fmt_head = "Fitness" + fmt_head_star * self.sol_size
        self.fmt_results = "{:7.2f}" + fmt_data_star * self.sol_size
        self.fmt_hbar = "-" * len(self.fmt_head)

    def print_empty_table(self):
        if self.n_top <= 50:
            for i in range(self.n_top + 6):
                print("")
        else:
            print("")

    def print_update(self):
        n_top = min(self.n_top, len(self.f))
        if n_top > 0:
            for i in range(6 + n_top):
                sys.stdout.write("\x1b[A")
            remain = self.elapsed / self.frac_done - self.elapsed
            if remain > 0:
                remain = time2human(remain)
            else:
                remain = "Done."
            print(
                "{percent:6.2f} %    Rate: {rate:>9,g}/s   Elapsed: {elapsed:<8}   ETA: {remain:<8}".format(
                    percent=min(1, self.frac_done) * 100,
                    rate=self.n_solved / self.elapsed,
                    elapsed=time2human(self.elapsed),
                    remain=remain,
                )
            )
            print(
                f"Gen: {self.gen:<6d} Best: {self.f[0]:>6.2f}        Average: {np.mean(self.f):>6.2f}   Worst: {self.f[-1]:>6.2f}"
            )
            print(self.fmt_hbar)
            print(self.fmt_head)
            print(self.fmt_hbar)
            for i, solution in enumerate(self.s[:n_top]):
                vals = list()
                vals.append(self.f[i])
                for j in range(self.sol_size):
                    index, offset = solution[j]
                    db_idx = self.db_idx[index]
                    dbindex = index - self.db_off[db_idx]
                    if self.db_n > 1:
                        vals.append(self.db_lab[db_idx])
                    vals.append(dbindex)
                    vals.append(np.log10(offset))
                print(self.fmt_results.format(*vals))
            print(self.fmt_hbar)
        else:
            sys.stdout.write("\x1b[A")
            sys.stdout.write("\x1b[A")
            print(
                "{percent:6.2f} %   Rate: {rate:>9,g}/s   Elapsed: {elapsed:<8}   ETA: {remain:<8}".format(
                    percent=self.frac_done * 100,
                    rate=self.n_solved / self.elapsed,
                    elapsed=time2human(self.elapsed),
                    remain=time2human(self.elapsed / self.frac_done - self.elapsed),
                )
            )

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
        if self.spread:
            mut_group = self.rng.integers(self.group_n, size=w.shape)
            mutants = self.group_off[mut_group] + self.rng.integers(
                self.group_num[mut_group]
            )
            mutants = self.group_index[mutants]
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
            w0 = w["offset"] * np.exp(
                self.rng.normal(scale=self.mut_offset_magnitude, size=w.shape)
            )
            w0 = w0 / (1 + w0)
            w["offset"] = np.choose(
                self.rng.uniform(size=w.shape) < self.mut_rate_offset,
                [
                    w["offset"],
                    w0,
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
        fc = get_fitness(
            self.trimmed_db,
            self.eval_data,
            self.exclude_index,
            c,
            fixed_offsets=self.fixed_offsets,
            ejecta=self.ejecta,
            cdf=self.cdf,
            limit_solution=self.limit_solution,
            limit_solver=self.limit_solver,
            local_search=self.local_search,
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

        # eliminate solutions that do not include pin'ned groups at correct numbers
        if self.pin_n > 0:
            igr = np.sort(self.group_idx[self.db_idx[o["index"]]], axis=-1)
            mask = np.full(igr.shape[0], True)
            for p, c in zip(self.pin, self.pin_count):
                mask &= np.sum(igr[:, :] == p, axis=1) == c
            o = o[mask]
            f = f[mask]
            if self.debug:
                n = np.sum(~mask)
                if n > 0:
                    self.logger.info(
                        f"pin: eliminating {n} candidates, {f.shape[0]} remaining."
                    )
            assert f.shape[0] > 0, "pin: no data left after candidate elimination"

        # If requested, eliminate solutions that do not span all groups.
        if self.spread:
            igr = np.sort(self.group_idx[self.db_idx[o["index"]]], axis=-1)
            for p in self.pin:
                igr = igr[igr != p].reshape(igr.shape[0], -1)
            if self.sol_size - self.pin_n >= self.pin_free_n:
                mask = np.all(igr[:, 1:] <= igr[:, :-1] + 1, axis=-1)
                mask &= igr[:, 0] == self.pin_free_min
                mask &= igr[:, -1] == self.pin_free_max
            else:
                mask = np.all(igr[:, 1:] > igr[:, :-1], axis=-1)
            o = o[mask]
            f = f[mask]
            if self.debug:
                n = np.sum(~mask)
                if n > 0:
                    self.logger.info(
                        f"spread: eliminating {n} candidates, {f.shape[0]} remaining."
                    )
            assert f.shape[0] > 0, "spread: no data left after candidate elimination"

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
        o = o[sel]
        f = f[sel]
        if self.debug:
            n = np.sum(~sel)
            if n > 0:
                self.logger.info(f"eliminating {n} duplicates, {f.shape[0]} remaining.")
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
            if self.debug:
                n = len(discard_index)
                if n > 0:
                    self.logger.info(
                        f"eliminating {n} solutions, {f.shape[0]} remaining."
                    )
            assert f.shape[0] > 0, "no data left after elimiation of solutions"
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

        j0 = 0
        if self.pin_n > 0:
            for p, c in zip(self.pin, self.pin_count):
                ncomb = int(comb(self.group_num[p], c))
                if ncomb > self.pop_size:
                    i0 = self.rng.integers(ncomb - self.pop_size)
                    idx = gen_slice(i0, i0 + self.pop_size, c)
                    ii = self.rng.permutation(self.pop_size)
                    idx = idx[ii]
                else:
                    idx = np.ndarray((self.pop_size, c), dtype=np.int64)
                    n = 0
                    while n < self.pop_size:
                        m = min(ncomb, self.pop_size - n)
                        ii = self.rng.permutation(m)
                        idx[n : n + m, :] = gen_slice(0, m, c)[ii]
                        n += m
                idx = idx + self.group_off[p]
                idx = self.group_index[idx]
                j1 = j0 + c
                s["index"][:, j0:j1] = idx
                j0 = j1

        if self.spread:
            assert np.all(
                self.group_num > (self.sol_size - self.group_n + 1)
            ), "at least one group has too few elements for requested solution size"
            k = self.sol_size - self.pin_n
            igr = self.rng.permuted(
                np.tile(
                    np.arange(self.pin_free_n, dtype=np.int64),
                    (self.pop_size, 1),
                ),
                axis=-1,
            )[:, :k]

            n = self.sol_size - self.pin_n - self.pin_free_n
            if n > 0:
                igr = np.concatenate(
                    (
                        igr,
                        self.rng.integers(
                            self.pin_free_n,
                            size=(
                                self.pop_size,
                                n,
                            ),
                        ),
                    ),
                    axis=-1,
                )

            igr = self.pin_free[igr]
            idx = self.group_off[igr] + self.rng.integers(self.group_num[igr])

            # replace duplicates
            while True:
                idx_sorted = np.sort(idx, axis=-1)
                ii = np.any(idx_sorted[:, 1:] == idx_sorted[:, :-1], axis=-1)
                if not np.any(ii):
                    break
                idx[ii] = self.group_off[igr[ii]] + self.rng.integers(
                    self.group_num[igr[ii]]
                )
            s["index"][:, j0:] = self.group_index[idx]
        else:
            # have few duplicates
            n = 0
            while n < self.pop_size:
                k = self.sol_size - self.pin_n
                m = min(self.db_size // k, self.pop_size - n)
                s["index"][n : n + m, j0:] = self.rng.permutation(self.db_size)[
                    : m * k
                ].reshape(-1, k)
                n += m
        s["offset"] = (
            self.rng.uniform(size=(self.pop_size, self.sol_size))
            * (0.9 / self.sol_size)
            + 1.0e-14
        )

        return s

    def plot_fitness(
        self,
        gen=False,
        fig=None,
        ax=None,
        show_copyright=True,
        figsize=(10, 6),
        dpi=102,
    ):
        """
        plot fitness as a function of time
        """

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
        if gen:
            ax.set_xlabel("Generations")
            times = np.arange(self.gen + 1)
        else:
            ax.set_xlabel("Time (s)")
            times = self.times

        ax.set_ylabel("Fitness (Error)")
        leg_starname(ax, self.star.name)
        if show_copyright:
            leg_copyright(ax)
        info = f"Generations: {self.gen:d}\nPopulation size: {self.pop_size:d}\nGene size: {self.sol_size:d}"
        leg_info(ax, info)

        ax.plot(times, self.history["worst"], label="worst", color="tab:red")
        ax.plot(times, self.history["mean"], label="mean", color="tab:blue")
        ax.plot(times, self.history["median"], label="median", color="tab:cyan")
        ax.plot(times, self.history["best"], label="best", color="tab:green")
        ax.set_yscale("log")

        leg = ax.legend(loc=(0.5, 0.1))
        leg.set_draggable(True)

        if fig is not None:
            fig.tight_layout()
