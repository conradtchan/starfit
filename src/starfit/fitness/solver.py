import numpy as np

from ._solver import fitness_


def check_offsets(offsets):
    if np.any(offsets > 1):
        ii = np.any(offsets > 1, axis=-1)
        i = np.where(ii)[0][0]
        raise AttributeError(f"invalid offset > 1 at index {i}: {offsets[i]}")


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
    ls=False,
):

    nsol = abu.shape[0]
    nstar = abu.shape[1]
    nel = abu.shape[2]
    ncov = covariance.shape[1]

    if cdf is None:
        cdf = True
    if cdf:
        icdf = 1
    else:
        icdf = 0

    check_offsets(offsets)

    ils = LS_MODES.get(ls, None)
    if ils is None:
        raise Exception(f"Unknown search option {ls=}")

    fitness, offset = fitness_(
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
    )

    check_offsets(offset)

    return fitness, offset


def testf():
    # obs = [-5.0, -5.0]
    # err = [0.3, -0.15]
    # abu = [
    #     [1.0e-05, 1.0e-05],
    #     [1.7e-05, 1.3e-05],
    # ]
    # nstar = 2
    # nel = 2

    # obs = [ -5.86418771e0, -6.19418771e0, -6.23418771e0, -9.26418771e0, -8.39418771e0, -9.48418771e0, -10.60418771e0, -12.36418771e0, -10.20418771e0 ]
    # err = [0.24e0, 0.3e0, 0.24e0, 0.08e0, 0.07e0, 0.09e0, 0.2e0, 0.2e0, 0.15e0]
    # abu = [
    # [ 9.23837878e-04, 2.86915484e-08, 9.86329569e-04, 7.22169530e-07, 1.44506086e-05, 2.98930995e-07, 3.79484481e-13, 3.49929099e-15, 9.71069640e-18],
    # [ 1.29517519e-03, 3.96547937e-07, 7.74776042e-03, 5.80834667e-06, 1.79736058e-04, 2.76109507e-06, 1.32311697e-05, 8.84902251e-08, 1.15314210e-05],
    # [ 1.06639960e-06, 1.63124275e-05, 3.58418401e-07, 8.92132488e-12, 1.22640917e-11, 2.84983024e-13, 8.85979959e-13, 5.86733432e-16, 1.37888246e-17],
    # ]
    # nstar = 3
    # nel = 9

    # y1 = []
    # y2 = []
    # y3 = []
    # x = []
    # for i in range(500):
    #     c = [0.1, i*0.01]
    #     f, f1, f2 = _solver.chisq(c, obs, err, abu, nstar, nel)
    #     x += [c[1]]
    #     y1 += [f]
    #     y2 += [f1[1]]
    #     y3 += [f2[1]]
    # plot(x,y1)
    # plot(x,y2)
    # plot(x,y3)

    # === interface has changed ===
    # c = np.ndarray((3,))
    # c[:] = 1.0e-7
    # cnew = _solver.solver.newton(c, obs, err, abu, nstar, nel)
    # print(cnew)
    # f, f1, f2 = _solver.solver.chisq(cnew, obs, err, abu, nstar, nel)
    # print(f)

    pass
