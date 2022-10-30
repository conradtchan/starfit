import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from uuid import uuid1

import numpy as np
from scipy.special import comb

from .autils.human import time2human
from .autils.utils import is_iterable
from .fit import get_solution
from .starfit import StarFit
from .utils import set_priority


class Multi(StarFit):
    """
    Find the best fit for multiple stars (complete search).
    The search results are printed live.

    partition:
        if number of data bases matches sol_size, take one from each DB

    sol_size:
        None:
            if partition is True, set sol_size to number of DBs

    TODO - allow multiple contributions form one partition
           (sol_size becomes vector with length of number of partitions)
    """

    def __init__(
        self,
        *args,
        n_top=10,  # How many solutions do we want to see?
        fixed_offsets=False,  # Fixed offsets based on ejecta mass?
        threads=None,
        nice=19,
        save=False,
        webfile=None,
        path=None,
        block_size=2**17,
        sol_size=None,
        partition=None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        if sol_size is None:
            if partition is None:
                partition = True
            if partition is True:
                sol_size = self.db_n
            elif partition is False:
                # make a million matches
                sol_size = max(1, int(6 / np.log10(self.db_size)))
                self.logger.info(f"setting {sol_size=}")
            elif is_iterable(partition):
                sol_size = len(partition)
                self.logger.info(f"setting {sol_size=}")
        if partition is None:
            partition = sol_size == self.db_n
        if sol_size == 1:
            partition = False
        if partition is True:
            partition = [[i] for i in range(self.db_n)]
        elif partition is False:
            partition = [[i for i in range(self.db_n)]]
        elif is_iterable(partition):
            if np.all([isinstance(i, int) for i in partition]):
                pnew = list()
                i = 0
                for p in partition:
                    assert p > 0, "invalid partition size {p=}"
                    pnew.append([j + i for j in range(p)])
                    i += p
                assert (
                    i == self.db_n
                ), f"invalid total number of DBs in partition n={i} vs. db_n={self.db_n}"
                partition = pnew
            else:
                pt = list()
                for p in partition:
                    assert is_iterable(p)
                    for pi in p:
                        assert isinstance(pi, int)
                        assert pi >= 0
                        assert pi < self.db_n
                        assert pi not in pt
                        pt.append(p)
                assert len(pt) == self.db_n
                assert len(set(pt)) == self.db_n
        else:
            raise AttributeError(f"[{self.__class__.__name__}] unknown {partition=}")

        partition_n = len(partition)
        partition_num = np.array([np.sum(self.db_num[p]) for p in partition])

        if is_iterable(sol_size):
            sol_sizes = np.array(sol_size, dtype=np.int64)
            sol_size = int(np.sum(sol_sizes))
        else:
            if partition_n == 1:
                sol_sizes = np.array([sol_size])
            else:
                sol_sizes = np.full(sol_size, 1, dtype=np.int64)

        if len(partition) != len(sol_sizes):
            raise AttributeError(f"require {len(partition)=} == {len(sol_sizes)=}")

        partition_comb = np.array(
            [int(comb(p, s)) for p, s in zip(partition_num, sol_sizes)]
        )

        partition_index = np.ndarray(self.db_size, dtype=np.int64)
        i0 = 0
        for pp in partition:
            for p in pp:
                n = self.db_num[p]
                i1 = i0 + n
                partition_index[i0:i1] = np.arange(n) + self.db_off[p]
                i0 = i1

        self.sol_size = sol_size
        self.sol_sizes = sol_sizes
        self.partition = partition
        self.partition_n = partition_n
        self.partition_num = partition_num
        self.partition_comb = partition_comb
        self.partition_index = partition_index

        self.fixed_offsets = fixed_offsets
        self.n_top = n_top

        self.n_combinations = np.product(self.partition_comb)

        if self.n_top is None:
            self.n_top = min(self.n_combinations, 20)

        self.block_size = min(block_size, self.n_combinations)

        self.logger.info(
            f"Searching {self.db_size:,d} models ({self.n_combinations:,d} combinations)"
        )

        if not self.silent:
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

            if not self.silent:
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
            initializer=set_priority,
            initargs=(nice,),
        )

        n_blocks = int(np.ceil(self.n_combinations / self.block_size))
        slice_range = np.arange(n_blocks + 1, dtype="int") * self.block_size
        slice_range[-1] = self.n_combinations

        futures = list()
        for i in range(n_blocks):
            futures.append(
                executor.submit(
                    get_solution,
                    gen_start=slice_range[i],
                    gen_end=slice_range[i + 1],
                    fixed=self.fixed_offsets,
                    eval_data=self.eval_data,
                    exclude_index=self.exclude_index,
                    trimmed_db=self.trimmed_db,
                    ejecta=self.ejecta,
                    sol_size=self.sol_size,
                    cdf=self.cdf,
                    num=self.partition_num,
                    size=self.sol_sizes,
                    com=self.partition_comb,
                    index=self.partition_index,
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
        if len(self.partition) > 1:
            print(f"Partitioning: {' x '.join(str(i) for i in self.partition_comb)}")
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
                        vals.append(db_idx)
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
            if self.fixed_offsets:
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
