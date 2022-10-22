"""Fitting stars"""

import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import getpid
from pathlib import Path
from uuid import uuid1

import numpy as np
import psutil
from scipy.special import comb

from .autils.human import time2human
from .autils.logged import Logged
from .solgen._solgen import gen_slice
from .starfit import Results, _fitness


class Single(Results, Logged):
    """Find the best fit for single stars (complete search)"""

    sol_size = 1

    def __init__(
        self,
        filename,
        db="znuc2012.S4.star.el.y.stardb.gz",
        *,
        combine=[[]],
        z_exclude=[],
        z_min=1,
        z_max=30,
        z_lolim=[],
        upper_lim=True,
        silent=False,
        cdf=True,
    ):
        self.silent = silent
        super().__init__()

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

        self.gene_size = 1
        # Initialize array
        stars = np.recarray(
            (self.db_size, 1), dtype=[("index", np.int64), ("offset", np.float64)]
        )
        stars.offset = 1.0e-4
        stars.index[:, 0] = np.arange(self.db_size)

        t_start = time.perf_counter()
        fitness = _fitness(
            trimmed_db=self.trimmed_db,
            eval_data=self.eval_data,
            z_exclude_index=self.exclude_index,
            sol=stars,
            cdf=cdf,
            ls=True,
        )
        self.logger.info("Calculation time:")
        self.logger.info(f" {time2human(time.perf_counter() - t_start)}")

        # Identify the best fit by finding the entry with the smallest error
        self.unsorted_fitness = fitness
        sort_index = np.argsort(fitness)
        self.sorted_fitness = fitness[sort_index]
        self.sorted_stars = stars[sort_index]

        self.logger.info("Best solution:")
        self.logger.info(" " + str(self.sorted_stars[0, 0]))
        self.logger.info("Fitness:")
        self.logger.info(" " + str(self.sorted_fitness[0]))

        self.close_logger(timing="Finished in")


class Multi(Results, Logged):
    """
    Find the best fit for multiple stars (complete search).
    The search results are printed live.

    partition:
        if number of data bases matches sol_size, take one from each DB

    sol_size:
        None:
            if partition is True, set sol_size to number of DBs
    """

    def __init__(
        self,
        filename,
        db="znuc2012.S4.star.el.y.stardb.gz",
        n_top=10,  # How many solutions do we want to see?
        fixed=False,  # Fixed offsets based on ejecta mass?
        threads=None,
        nice=19,
        cdf=True,
        combine=[[]],
        z_min=1,
        z_max=30,
        z_exclude=[],
        z_lolim=[],
        upper_lim=True,
        silent=False,
        save=False,
        webfile=None,
        path=None,
        block_size=2**17,
        sol_size=None,
        partition=None,
    ):
        self.silent = silent
        super().__init__()

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
            if partition is None:
                partition = True
            if partition:
                sol_size = self.db_n
            else:
                # make a million matches
                sol_size = max(1, int(6 / np.log10(self.db_size)))
                self.logger.info(f"setting {sol_size=}")
        if partition is None:
            partiton = sol_size == self.db_n
        if sol_size == 1:
            partition = False
        if partition and sol_size != self.db_n:
            raise AttributeError(
                f"currently {partiton=} requires {sol_size=} == len(DB)={self.db_n}"
            )

        self.sol_size = sol_size
        self.partition = partition

        self.fixed = fixed
        self.n_top = n_top

        if self.partition:
            self.n_combinations = int(np.product(self.db_num))
        else:
            self.n_combinations = int(comb(self.db_size, self.sol_size))

        if self.n_top is None:
            self.n_top = self.n_combinations

        self.block_size = min(block_size, self.n_combinations)
        self.cdf = cdf

        self.logger.info(
            f"Searching {self.db_size:,d} models ({self.n_combinations:,d} combinations)"
        )

        if not silent:
            self.print_header()

        self.working_arrays()

        time_start = time.perf_counter()
        futures = self.init_futures(threads=threads, nice=nice)

        if not self.silent:
            self.print_empty_table()

        elapsed = 0
        self.n_solved = 0

        for n, future in enumerate(as_completed(futures)):
            (
                local_stars,
                local_fitness,
                n_local_solved,
            ) = future.result()
            return_len = len(local_stars)

            self.top_stars[2 * self.n_top - return_len :, :] = local_stars
            self.top_fitness[2 * self.n_top - return_len :] = local_fitness

            sort = np.argsort(self.top_fitness)
            self.top_stars = self.top_stars[sort, :]
            self.top_fitness = self.top_fitness[sort]

            self.n_solved += n_local_solved
            self.frac_done = self.n_solved / self.n_combinations

            elapsed = time.perf_counter() - time_start
            try:
                self.cycle_time = (
                    0.01 * (elapsed - self.elapsed) + 0.99 * self.cycle_time
                )
            except:
                self.cycle_time = 1.0e-10

            self.elapsed = elapsed

            if not silent:
                self.print_update()

        self.sorted_stars = self.top_stars[: self.n_top]
        self.sorted_fitness = self.top_fitness[: self.n_top]

        if save:
            if webfile:
                filename = webfile
                path = Path("/tmp")
            else:
                if save is True:
                    filename = uuid1().hex + ".txt"
                else:
                    filename = save
                if path is None:
                    path = "."
                path = Path(path).expanduser().resolve()
            file = path / filename
            self.save_file(file)

        self.close_logger(timing="Finished in")

    def working_arrays(self):
        sorting_size = 2 * self.n_top
        self.top_stars = np.ndarray(
            (sorting_size, self.sol_size),
            dtype=[("index", np.int64), ("offset", np.float64)],
        )
        self.top_fitness = np.ndarray((sorting_size,), dtype=np.float64)
        self.top_fitness[:] = np.inf

    def init_futures(self, threads=None, nice=19):
        if not self.silent:
            print("Initializing...")

        if threads is None:
            threads = multiprocessing.cpu_count()
        executor = ProcessPoolExecutor(
            max_workers=threads,
            initializer=_set_priority,
            initargs=(nice,),
        )

        n_blocks = int(np.ceil(self.n_combinations / self.block_size))
        slice_range = np.arange(n_blocks + 1, dtype="int") * self.block_size
        slice_range[-1] = self.n_combinations

        if self.partition:
            num = self.db_num
        else:
            num = None

        futures = list()
        for i in range(n_blocks):
            futures.append(
                executor.submit(
                    _solve,
                    gen_start=slice_range[i],
                    gen_end=slice_range[i + 1],
                    fixed=self.fixed,
                    eval_data=self.eval_data,
                    exclude_index=self.exclude_index,
                    trimmed_db=self.trimmed_db,
                    ejecta=self.ejecta,
                    sol_size=self.sol_size,
                    cdf=self.cdf,
                    num=num,
                    return_size=self.n_top,
                )
            )
            if not self.silent:
                sys.stdout.write("\x1b[A")
                print(f"{(i+1)/n_blocks:>6.2f} % initialized.")

        if not self.silent:
            sys.stdout.write("\x1b[A")
            print("Waiting for first results...")
        return futures

    def print_header(self):
        print("\n")
        print("=================================================================")
        print(
            f"Matching {self.n_combinations:,d} combinations of {self.sol_size} stars"
        )
        if self.partition:
            print(f"Partitioning: {' x '.join(str(i) for i in self.db_num)}")
        if self.db_n == 1:
            print(f"Database: {self.db[0].name}")
        else:
            print("Databases:")
            # for i, db in enumerate(self.db):
            #     print(f"{i+1:>6d}: {db.name}")
            self.print_db(ind=4)
        print(f"Star: {self.star.name}")
        print(
            f"Using {self.fit_size:d} elements ({np.sum(self.eval_data.error < 0):d} limits)"
        )
        if self.fixed:
            print("Offsets are ejecta mass")
        else:
            print("Offsets solved using Newton-Raphson")
        print("=================================================================")

        fmt_head_star = " Index  Offset"
        fmt_data_star = "{:>6d}  {:6.2f}"
        fmt_pad = " " * 5
        if self.db_n > 1:
            fmt_head_star = "DB " + fmt_head_star
            fmt_data_star = "{:>2d} " + fmt_data_star
        fmt_head_star = fmt_pad + fmt_head_star
        fmt_data_star = fmt_pad + fmt_data_star
        self.fmt_head = "Fitness" + fmt_head_star * self.sol_size
        self.fmt_results = "{:7.2f}" + fmt_data_star * self.sol_size
        self.fmt_hbar = "-" * len(self.fmt_head)

    def print_empty_table(self):
        sys.stdout.write("\x1b[A")
        if self.n_top <= 50:
            for i in range(self.n_top + 6):
                print("")
        else:
            print("")

    def print_update(self):
        if self.n_top <= 20:
            for i in range(6 + self.n_top):
                sys.stdout.write("\x1b[A")
            print(
                "{percent:6.2f} %    Rate: {rate:>9,g}/s   Elapsed: {elapsed:<8}   ETA: {remain:<8}".format(
                    percent=self.frac_done * 100,
                    rate=self.n_solved / self.elapsed,
                    elapsed=time2human(self.elapsed),
                    remain=time2human(self.elapsed / self.frac_done - self.elapsed),
                )
            )
            print("")
            print(self.fmt_hbar)
            print(self.fmt_head)
            print(self.fmt_hbar)
            for i, solution in enumerate(self.top_stars[: self.n_top]):
                vals = list()
                vals.append(self.top_fitness[i])
                for j in range(self.sol_size):
                    index, offset = solution[j]
                    db_idx = self.db_idx[index]
                    dbindex = index - self.db_off[db_idx]
                    if self.db_n > 1:
                        vals.append(db_idx + 1)
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

    def save_file(self, file_name):
        file = Path(file_name).expanduser().resolve()
        with open(file, mode="wt") as f:
            text = self.format(n=len(self.top_stars), wide=99)
            lines = text.splitlines()
            totlen = np.max([len(line) for line in lines])
            barlen = totlen + 3
            hbar = "-" * barlen + "\n"
            xbar = "=" * barlen + "\n"
            f.write(xbar)

            f.write(
                f"Matching {self.n_combinations:,d} combinations of {self.sol_size:d} stars\n"
            )
            if self.db_n == 1:
                f.write(f"Database: {self.db[0].name}\n")
            else:
                f.write("Databases:\n")
                for l in self.format_db(ind=4).splitlines():
                    f.write(l + "\n")
            f.write(f"Star: {self.star.name}\n")
            f.write(
                f"Using {self.fit_size:d} elements ({np.sum(self.eval_data.error < 0):d} upper limits)\n"
            )
            if self.fixed:
                f.write("Offsets are according to ejecta mass\n")
            else:
                f.write("Offsets solved using UOBYQA\n")

            f.write(xbar)
            f.write(hbar)
            for l in lines[:2]:
                f.write(l + "\n")
            f.write(hbar)
            for l in lines[2:]:
                f.write(l + "\n")
            f.write(hbar)

            # write db comments / info
            for l in self.format_comments().splitlines():
                f.write(l + "\n")


def _set_priority(value: int):
    p = psutil.Process(getpid())
    p.nice(value)


def gen_map(gen_start, gen_end, num):
    off = np.cumsum(num)
    off[1:] = off[0:-1]
    off[0] = 0
    idx = np.arange(gen_start, gen_end, dtype=np.int64)
    idx = np.array(np.unravel_index(idx, num)).T
    return idx + off


def _solve(
    gen_start,
    gen_end,
    exclude_index,
    trimmed_db,
    fixed=False,
    ejecta=None,
    return_size=None,
    sol_size=None,
    num=None,
    **kwargs,
):
    """
    Local search solve for a set of stars. Used for splitting into different
    processes. This needs to be a separate function that can be pickled.
    """

    solutions = np.ndarray(
        (gen_end - gen_start, sol_size), dtype=[("index", "int"), ("offset", "f8")]
    )

    if num is None:
        solutions["index"][:] = gen_slice(gen_start, gen_end, sol_size)
    else:
        solutions["index"][:] = gen_map(gen_start, gen_end, num)

    if fixed:
        # for i in range(2):
        #     solutions["offset"][:, i] = ejecta[solutions["index"][:, i]]
        solutions["offset"][:, :] = ejecta[solutions["index"][:, :]]
    else:
        solutions[:, :]["offset"] = 1.0e-5

    fitness = _fitness(
        sol=solutions,
        trimmed_db=trimmed_db,
        fixed_offsets=fixed,
        ejecta=ejecta,
        ls=not fixed,
        z_exclude_index=exclude_index,
        **kwargs,
    )

    sort = np.argsort(fitness)
    solutions = solutions[sort, :]
    fitness = fitness[sort]

    n_local_solved = len(fitness)

    if return_size is not None:
        solutions = solutions[:return_size, :]
        fitness = fitness[:return_size]

    return solutions, fitness, n_local_solved


class Double(Multi):
    sol_size = 2

    def __init__(self, *args, **kwargs):
        sol_size = kwargs.setdefault("sol_size", self.sol_size)
        assert sol_size == 2, f"require sol_size == 2, provided: {sol_size=}"
        super().__init__(*args, **kwargs)


class Direct(Results, Logged):
    """Find the best fit for stars you specify"""

    def __init__(
        self,
        filename,
        db="znuc2012.S4.star.el.y.stardb.gz",
        stars=[[114, 3604]],  # A list of the indices of the stars
        offsets=None,
        z_max=30,
        fig=None,
        silent=False,
        combine=[[]],
        z_exclude=[],
        z_lolim=[],
        upper_lim=True,
        fixed=False,
        cdf=True,
    ):
        self.silent = silent
        Results.__init__(self, "Direct")

        self._setup(
            filename,
            db,
            combine,
            z_exclude,
            z_max,
            upper_lim,
            z_lolim,
        )

        if stars is not None:
            stars = np.array(stars)
            self.cdf = cdf
            sol, fitness = self.run(stars, offsets, fixed)

            sort = np.argsort(fitness)
            self.sorted_fitness = fitness[sort]
            self.sorted_stars = sol[sort, :]

            self.logger.info("Stars fitted:")
            for i in range(stars.shape[0]):
                self.logger.info(f"Solution {i}:")
                for j in range(stars.shape[1]):
                    self.logger.info(f" {stars[i, j]:05}:    {sol[i, j]['offset']:5e}")

                self.logger.info(" Fitness:")
                self.logger.info(f"  {fitness[i]:7.3f}")

            self.close_logger(timing="Finished in")
