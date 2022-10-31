import numpy as np

import starfit


def test():
    """Test a full binary search"""
    TEST_STAR = "HE1327-2326.dat"
    TEST_DB = "he2sn.HW02.star.el.y.stardb.gz"
    EXPECTED_RESULT = 32.885971124847515

    test_result = starfit.Multi(
        TEST_STAR, TEST_DB, sol_size=2, group=False, silent=True, z_max=30
    ).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
