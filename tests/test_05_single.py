import numpy as np

import starfit


def test():
    """Test a single star fits, data format 5"""
    TEST_STAR = "CI-Chondrite-Lodders2019_f5.dat"
    EXPECTED_RESULT = 16372.053997257417

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
