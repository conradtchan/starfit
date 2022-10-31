import numpy as np

import starfit


def test():
    """Test a partitioned multi-star multi-DB search"""
    TEST_STAR = "HE1327-2326.dat"
    TEST_DB = (
        "he2sn.HW02.star.el.y.stardb.gz",
        "znuc_lc12.el.y.stardb.xz",
        "he2sn.hennuc.star.el.y.stardb.xz",
    )
    EXPECTED_RESULT = 15.337664935301007

    test_result = starfit.Multi(
        TEST_STAR, TEST_DB, silent=True, z_max="Zn", group=True
    ).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
