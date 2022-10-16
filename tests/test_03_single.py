import numpy as np

import starfit


def test():
    """Test a single star fits, data format 3"""
    TEST_STAR = "BD_80_245.dat"
    EXPECTED_RESULT = 1.12390641724435

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
