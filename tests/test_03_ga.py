from numpy import pi as PI

import starfit


def test():
    """Test a single star fit"""
    TEST_STAR = "HE1327-2326.dat"
    FITNESS_THRESHOLD = PI

    test_result = starfit.Ga(
        TEST_STAR,
        silent=True,
        time_limit=1,
        sol_size=3,
        z_max=30,
    ).sorted_fitness[0]

    assert test_result < FITNESS_THRESHOLD
