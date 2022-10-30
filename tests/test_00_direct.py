import numpy as np

import starfit


def test():
    """Test a direct single star match"""
    TEST_STAR = "HE1327-2326.dat"
    TEST_DB = "he2sn.HW02.star.el.y.stardb.gz"
    EXPECTED_RESULT = 34.188915103255994
    TEST_MODELS = [[0, 1, 2]]

    test_result = starfit.Direct(
        TEST_STAR, TEST_DB, stars=TEST_MODELS, silent=True, z_max=30
    ).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
