"""Fitting stars"""

import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import psutil

from .autils.human import time2human
from .autils.logged import Logged
from .solgen._solgen import comb, gen_slice
from .starfit import Results


class Single(Results, Logged):
    """Find the best fit for single stars (complete search)"""

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
        super().__init__("Single")

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
        db_size = self.trimmed_db.shape[1]
        # Initialize array
        stars = np.recarray((db_size, 1), dtype=[("index", "int"), ("offset", "f8")])
        stars.offset = 1.0e-4
        stars.index[:, 0] = np.arange(db_size)

        t_start = time.perf_counter()
        fitness = self._fitness(
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


class Double(Results, Logged):
    """
    Find the best fit for 2 stars (complete search).
    The search results are printed live.

    TODO - multiple stars
    TODO - 'partition' to select one from each db
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
        block_size=2**17,
    ):
        self.silent = silent
        Results.__init__(self, "Double")

        self._setup(
            filename=filename,
            dbname=db,
            combine=combine,
            z_exclude=z_exclude,
            z_min=z_min,
            z_max=z_max,
            upper_lim=upper_lim,
            z_lolim=z_lolim,
        )

        self.setup_double(fixed, n_top, block_size, cdf)

        self.logger.info(
            f"Searching {self.db_size:_d} models ({self.n_solves:_d} combinations)"
        )

        if not silent:
            self.print_header()

        self.working_arrays()

        time_start = time.perf_counter()
        futures = self.init_futures(threads=threads, nice=nice)

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
            self.frac_done = self.n_solved / self.n_solves

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
                directory = "/tmp"
            else:
                if fixed:
                    string = " (fixed)"
                else:
                    string = ""
                filename = "{:s}({:s}, {:s}) {:d}.{:d}{:s}".format(
                    self.__class__.__name__,
                    self.star.name,
                    "-".join([db.name for db in self.db]),
                    self.n_top,
                    self.db_size,
                    string,
                )
                path = os.path.dirname(__file__)
                directory = os.path.join(path, "../results")
            file = os.path.join(directory, filename)
            self.save_file(file)

        self.close_logger(timing="Finished in")

    def setup_double(self, fixed, n_top, block_size, cdf):
        self.fixed = fixed
        self.n_top = n_top

        self.db_size = self.trimmed_db.shape[1]

        self.n_solves = self.db_size * (self.db_size - 1) // 2
        if self.n_top is None:
            self.n_top = self.n_solves

        self.block_size = min(block_size, self.n_solves)
        self.cdf = cdf

    def working_arrays(self):
        sorting_size = 2 * self.n_top
        self.top_stars = np.ndarray(
            (sorting_size, 2), dtype=[("index", "int"), ("offset", "f8")]
        )
        self.top_fitness = np.ndarray((sorting_size,), dtype="f8")
        self.top_fitness[:] = np.inf

    def init_futures(self, threads=None, nice=19):
        if threads is None:
            threads = multiprocessing.cpu_count()
        executor = ProcessPoolExecutor(
            max_workers=threads,
            initializer=_set_priority,
            initargs=(nice,),
        )

        n_combinations = int(comb(self.db_size, self.sol_size))
        n_blocks = int(np.ceil(n_combinations / self.block_size))
        slice_range = np.arange(n_blocks + 1, dtype="int") * self.block_size
        slice_range[-1] = n_combinations

        futures = [
            executor.submit(
                _solve,
                gen_start=slice_range[i],
                gen_end=slice_range[i + 1],
                fixed=self.fixed,
                eval_data=self.eval_data,
                exclude_index=self.exclude_index,
                trimmed_db=self.trimmed_db,
                ejecta=self.ejecta,
                cdf=self.cdf,
                return_size=self.n_top,
            )
            for i in range(n_blocks)
        ]

        if not self.silent:
            sys.stdout.write("\x1b[A")

        return futures

    def print_header(self):
        # TODO - update for db
        # TODO - update for multistar
        print("\n")
        print("=================================================================")
        print(
            "Matching {count:_d} combinations of 2 stars\nDatabase: {dbname}\nStar: {starname}\nUsing {nel:d} elements ({nup:d} upper limits)".format(
                count=self.db_size * (self.db_size - 1) // 2,
                dbname=", ".join([db.name for db in self.db]),
                starname=self.star.name,
                nel=self.trimmed_db.shape[0],
                nup=np.sum(self.eval_data.error < 0),
            )
        )
        if self.fixed:
            print("Offsets are ejecta mass")
        else:
            print("Offsets solved using Newton-Raphson")
        print("=================================================================")

        print("Initializing...")
        if self.n_top <= 50:
            for i in range(self.n_top + 7):
                print("")
        else:
            print("")

    def print_update(self):
        # TODO - update for DB
        # TODO - update for multistar
        if self.n_top <= 20:
            for i in range(7 + self.n_top):
                sys.stdout.write("\x1b[A")
            print(
                "{percent:2.2f} %     Rate: {rate:<6.0f}/s     Elapsed: {elapsed:<8}     ETA: {remain:<8}".format(
                    percent=self.frac_done * 100,
                    rate=self.block_size / self.cycle_time,
                    elapsed=time2human(self.elapsed),
                    remain=time2human(self.elapsed / self.frac_done - self.elapsed),
                )
            )
            print("")
            print("-----------------------------------------------------------------")
            print("Index    Offset            Index    Offset             Fitness")
            print("-----------------------------------------------------------------")
            for i, solution in enumerate(self.top_stars[: self.n_top]):
                print(
                    "{index1:<5}    {offset1:e}      {index2:<5}    {offset2:e}      {fitness:5.2f}".format(
                        index1=solution["index"][0],
                        index2=solution["index"][1],
                        offset1=solution["offset"][0],
                        offset2=solution["offset"][1],
                        fitness=self.top_fitness[i],
                    )
                )
            print("-----------------------------------------------------------------")
        else:
            sys.stdout.write("\x1b[A")
            sys.stdout.write("\x1b[A")
            print(
                "{percent:2.2f} %     Rate: {rate:<6.0f}/s     Elapsed: {elapsed:<8}     ETA: {remain:<8}".format(
                    percent=self.frac_done * 100,
                    rate=self.n_solved / self.elapsed,
                    elapsed=time2human(self.elapsed),
                    remain=time2human(self.elapsed / self.frac_done - self.elapsed),
                )
            )

    def save_file(self, file_name):
        # TODO - update for multistar
        with open(file_name, mode="w") as f:
            # columns = [
            #     (
            #         "Index",
            #         "{index1:<5}",
            #         5,
            #     ),
            #     (
            #         "Mass",
            #         "{mass1:5.1f}",
            #         5,
            #     ),
            #     (
            #         "Eexp",
            #         "{energy1:4.1f}",
            #         4,
            #     ),
            #     (
            #         "Mixing",
            #         "{mixing1:8.6f}",
            #         8,
            #     ),
            #     (
            #         "Offset",
            #         "{offset1:12e}",
            #         12,
            #     ),
            #     (
            #         "|",
            #         "|",
            #         1,
            #     ),
            #     (
            #         "Index",
            #         "{index2:<5}",
            #         5,
            #     ),
            #     (
            #         "Mass",
            #         "{mass2:5.1f}",
            #         5,
            #     ),
            #     (
            #         "Eexp",
            #         "{energy2:4.1f}",
            #         4,
            #     ),
            #     (
            #         "Mixing",
            #         "{mixing2:8.6f}",
            #         8,
            #     ),
            #     (
            #         "Offset",
            #         "{offset2:12e}",
            #         12,
            #     ),
            #     (
            #         "|",
            #         "|",
            #         1,
            #     ),
            #     (
            #         "Fitness",
            #         "{fitness:7.5f}",
            #         7,
            #     ),
            # ]

            # totlen = sum([c[2] for c in columns]) + len(columns) - 1
            # barlen = totlen + 3

            # f.write("=" * barlen + "\n")

            text = self.format(n=len(self.top_stars), wide=99)
            lines = text.splitlines()
            totlen = np.max([len(line) for line in lines])
            barlen = totlen + 3
            f.write("=" * barlen + "\n")

            f.write(
                "Matching {count:d}^2 combinations of 2 stars\nDatabase: {dbname}\nStar: {starname}\nUsing {nel:d} elements ({nup:d} upper limits)\n".format(
                    count=self.db_size,
                    dbname=", ".join([db.name for db in self.db]),
                    starname=self.star.name,
                    nel=self.trimmed_db.shape[0],
                    nup=np.sum(self.eval_data.error < 0),
                )
            )
            if self.fixed:
                f.write("Offsets are according to ejecta mass\n")
            else:
                f.write("Offsets solved using UOBYQA\n")

            f.write("=" * barlen + "\n")
            f.write("-" * barlen + "\n")
            for l in lines[:2]:
                f.write(l)
            f.write("-" * barlen + "\n")
            for l in lines[2:]:
                f.write(l)
            f.write("-" * barlen + "\n")

            # f.write("=" * barlen + "\n")
            # f.write("-" * barlen + "\n")
            # f.write(" ".join("{:<{}s}".format(c[0], c[2]) for c in columns) + "\n")
            # f.write("-" * barlen + "\n")
            # for i, solution in enumerate(self.top_stars):
            #     index1 = solution["index"][0]
            #     index2 = solution["index"][1]
            #     f.write(
            #         (" ".join([c[1] for c in columns]) + "\n").format(
            #             index1=index1,
            #             index2=index2,
            #             offset1=solution["offset"][0],
            #             offset2=solution["offset"][1],
            #             mass1=self.db.fielddata[index1]["mass"],
            #             mass2=self.db.fielddata[index2]["mass"],
            #             energy1=self.db.fielddata[index1]["energy"],
            #             energy2=self.db.fielddata[index2]["energy"],
            #             mixing1=self.db.fielddata[index1]["mixing"],
            #             mixing2=self.db.fielddata[index2]["mixing"],
            #             fitness=self.top_fitness[i],
            #         )
            #     )
            # f.write("-" * barlen + "\n")


def _set_priority(value: int):
    p = psutil.Process(os.getpid())
    p.nice(value)


def _solve(
    gen_start,
    gen_end,
    exclude_index,
    trimmed_db,
    fixed=False,
    ejecta=None,
    return_size=None,
    sol_size=2,
    **kwargs,
):
    """
    Local search solve for a set of stars. Used for splitting into different
    processes. This needs to be a separate function that can be pickled.
    """

    solutions = np.ndarray(
        (gen_end - gen_start, sol_size), dtype=[("index", "int"), ("offset", "f8")]
    )

    solutions["index"][:] = gen_slice(gen_start, gen_end, sol_size)

    if fixed:
        # for i in range(2):
        #     solutions["offset"][:, i] = ejecta[solutions["index"][:, i]]
        solutions["offset"][:, :] = ejecta[solutions["index"][:, :]]
    else:
        solutions[:, :]["offset"] = 1.0e-5

    fitness = Results._fitness(
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
