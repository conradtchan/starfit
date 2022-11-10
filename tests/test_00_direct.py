import numpy as np
import pytest

import starfit

testdata = [
    [[[0]], False, 392.91323408409573],
    [[[0]], True, 392.9132340840956],
    [[[0, 1, 2]], False, 34.18889930088918],
    [[[0, 1, 2]], True, 34.188915103255994],
]


@pytest.mark.parametrize("stars, cdf, expected", testdata)
def test(stars, cdf, expected):
    """Test a direct triple star match"""
    TEST_STAR = "HE1327-2326.dat"
    TEST_DB = "he2sn.HW02.star.el.y.stardb.gz"

    test_result = starfit.Direct(
        TEST_STAR, TEST_DB, stars=stars, silent=True, z_max=30, cdf=cdf
    ).sorted_fitness[0]

    assert np.isclose(test_result, expected)
