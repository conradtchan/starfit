"""Fitting stars"""

import numpy as np

from .fitness import solver
from .solgen._solgen import gen_slice


def get_fitness(
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

    # TODO: for multiple calls, such as GA, stardb (and star) data
    # should be passed only once.

    if fixed_offsets:
        offset = ejecta[sol["index"]]
        ls = False
    else:
        offset = sol["offset"]

    eval_data = eval_data[~z_exclude_index]
    abu = np.transpose(trimmed_db[~z_exclude_index][:, sol["index"]], (1, 2, 0))

    fitness, offsets = solver.fitness(
        eval_data.abundance,
        eval_data.error,
        eval_data.corr,
        abu,
        offset,
        cdf=cdf,
        ls=ls,
    )

    sol["offset"] = offsets
    fitness /= eval_data.error.shape[0] - 1
    return fitness


def gen_map(gen_start, gen_end, num):
    off = np.cumsum(num)
    off[1:] = off[0:-1]
    off[0] = 0
    idx = np.arange(gen_start, gen_end, dtype=np.int64)
    idx = np.array(np.unravel_index(idx, num)).T
    return idx + off


def get_solution(
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

    fitness = get_fitness(
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
