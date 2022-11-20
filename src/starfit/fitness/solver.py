import numpy as np

from ._solver import fitness_, fitness_m_


def check_offsets(offsets, flags, tag):
    if flags & FLAGS_LIMIT_SOLUTION:
        if np.any((s := np.sum(offsets, axis=-1)) > 1):
            i = np.where(s > 1)[0]
            raise AttributeError(
                f"[fitness][{tag}] invalid sum of offsets for solution(s) {i} > 1: {s[i]}"
            )
    if flags & FLAGS_LIMITED_SOLVER:
        if np.any(offsets > 1):
            ii = np.any(offsets > 1, axis=-1)
            i = np.where(ii)[0][0]
            raise AttributeError(
                f"[fitness][{tag}] invalid offset > 1 at index {i}: {offsets[i]}"
            )


ALMOST_ONE = 1.0e0 - 1.0e-14

# from f90 module
FLAGS_LIMIT_SOLUTION_BIT = 0
FLAGS_LIMITED_SOLVER_BIT = 1

# flags as addable integers
FLAGS_LIMIT_SOLUTION = 2**FLAGS_LIMIT_SOLUTION_BIT
FLAGS_LIMITED_SOLVER = 2**FLAGS_LIMITED_SOLVER_BIT

# default flags
FLAGS_DEFAULT = FLAGS_LIMIT_SOLUTION + FLAGS_LIMITED_SOLVER

# local search options
LS_CHI_ONLY = -1
LS_OFF = 0
LS_ON = 1
LS_PSOLVE = 2

LS_MODES = {
    None: LS_CHI_ONLY,
    False: LS_OFF,
    True: LS_ON,
    Ellipsis: LS_PSOLVE,
}


def fitness(
    observed,
    error,
    detection,
    covariance,
    abu,
    offsets,
    cdf=None,
    local_search=True,
    limit_solution=None,
    limit_solver=None,
    return_matrix=False,
):

    nsol = abu.shape[0]
    nstar = abu.shape[1]
    nel = abu.shape[2]
    ncov = covariance.shape[1]

    if limit_solution is None:
        limit_solution = True
    if limit_solver is None:
        limit_solver = True

    flags = 0
    if limit_solution:
        flags += FLAGS_LIMIT_SOLUTION
    if limit_solver:
        flags += FLAGS_LIMITED_SOLVER

    if cdf is None:
        cdf = True
    if cdf:
        icdf = 1
    else:
        icdf = 0

    try:
        check_offsets(offsets, flags, "INPUT")
    except:
        if flags & FLAGS_LIMIT_SOLUTION:
            s = np.sum(offsets, axis=-1)
            ii = s > ALMOST_ONE
            offsets[ii, :] = offsets[ii, :] / s[ii, np.newaxis]
        else:
            raise

    ils = LS_MODES.get(local_search, None)
    if ils is None:
        raise Exception(f"Unknown search option {local_search=}")

    if return_matrix:
        fitness_x = fitness_m_
    else:
        fitness_x = fitness_

    fitness, offset = fitness_x(
        c=offsets,
        obs=observed,
        err=error,
        det=detection,
        cov=covariance,
        abu=abu,
        nel=nel,
        ncov=ncov,
        nstar=nstar,
        nsol=nsol,
        ls=ils,
        icdf=icdf,
        flags=flags,
    )

    check_offsets(offset, flags, "RESULT")

    return fitness, offset
