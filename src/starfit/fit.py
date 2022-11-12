"""Fitting stars"""

import numpy as np
from scipy.special import comb

from .fitness import solver
from .solgen._solgen import gen_slice


def get_fitness(
    trimmed_db,
    eval_data,
    z_exclude_index,
    sol,
    ejecta=[],
    fixed_offsets=False,
    local_search=True,
    **kwargs,
):
    """
    Evaluate the fitness of a set of solutions.
    If abundance is directly given in the case of smart GA, use that.
    """

    if fixed_offsets:
        offset = ejecta[sol["index"]]
        if local_search in (
            True,
            Ellipsis,
        ):
            local_search = False
    else:
        offset = sol["offset"]

    eval_data = eval_data[~z_exclude_index]
    abu = np.transpose(trimmed_db[~z_exclude_index][:, sol["index"]], (1, 2, 0))

    fitness, offsets = solver.fitness(
        eval_data.abundance,
        eval_data.error,
        eval_data.detection,
        eval_data.covariance,
        abu,
        offset,
        local_search=local_search,
        **kwargs,
    )

    sol["offset"] = offsets
    fitness /= eval_data.error.shape[0] - 1
    return fitness


def gen_map(gen_start, gen_end, num, size=None, com=None):
    if size is None:
        size = np.full(len(num), 1, dtype=np.int64)
    if com is None:
        com = np.array(
            [int(comb(n, s)) for n, s in zip(num, size)],
            dtype=np.int64,
        )
    idx = np.arange(gen_start, gen_end, dtype=np.int64)
    idx = np.array(np.unravel_index(idx, com)).T

    off = np.cumsum(num)
    off[1:] = off[0:-1]
    off[0] = 0

    off_ = list()
    for o, s in zip(off, size):
        off_.extend([o] * s)
    off = np.array(off_)

    size_off = np.array([0] + np.cumsum(size).tolist())
    ind = np.ndarray((gen_end - gen_start, sum(size)), dtype=np.int64)
    for i, s in enumerate(size):
        if s == 1:
            ind[:, size_off[i]] = idx[:, i]
            continue
        ii = np.where(idx[:, i] == 0)
        if len(ii) == 0:
            ix = gen_slice(idx[0, i], idx[-1, i] + 1, s)
            ind[:, size_off[i] : size_off[i + 1]] = ix[idx[:, i] - idx[0, i], :]
        elif len(ii) > 1 or (len(ii) == 1 and idx[0, i] <= idx[-1, i] + 1):
            ix = gen_slice(0, com[i], s)
            ind[:, size_off[i] : size_off[i + 1]] = ix[idx[:, i], :]
        else:
            i = ii[0][0]
            ix = gen_slice(idx[0, i], comb[i], s)
            ind[:i, size_off[i] : size_off[i + 1]] = ix[idx[:i, i] - idx[0, i], :]
            ix = gen_slice(0, idx[-1, i], s)
            ind[i:, size_off[i] : size_off[i + 1]] = ix[idx[i:, i], :]

    return ind + off


def get_solution(
    gen_start,
    gen_end,
    exclude_index,
    trimmed_db,
    fixed_offsets=False,
    ejecta=None,
    return_size=None,
    sol_size=None,
    num=None,
    size=None,
    com=None,
    index=None,
    **kwargs,
):
    """
    Local search solve for a set of stars. Used for splitting into different
    processes. This needs to be a separate function that can be pickled.
    """

    solutions = np.ndarray(
        (gen_end - gen_start, sol_size),
        dtype=[("index", np.int64), ("offset", np.float64)],
    )

    solutions["index"][:] = index[gen_map(gen_start, gen_end, num, size, com)]

    if fixed_offsets:
        solutions["offset"][:, :] = ejecta[solutions["index"][:, :]]
    else:
        solutions[:, :]["offset"] = 1.0e-5

    fitness = get_fitness(
        sol=solutions,
        trimmed_db=trimmed_db,
        fixed_offsets=fixed_offsets,
        ejecta=ejecta,
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
