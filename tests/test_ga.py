import starfit


def test():
    """Test a single star fit"""
    TEST_STAR = "HE1327-2326.dat"
    FITNESS_THRESHOLD = 3.0

    test_result = starfit.Ga(
        TEST_STAR, silent=True, time_limit=20, sol_size=3
    ).sorted_fitness[0]

    assert test_result < FITNESS_THRESHOLD
