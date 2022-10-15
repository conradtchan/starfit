import numpy as np

import starfit


def test():
    """Test a single star fit"""
    TEST_STAR = "HE1327-2326.dat"
    EXPECTED_RESULT = 3.292526847047224

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
