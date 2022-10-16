import numpy as np

import starfit


def test():
    """Test a single star fits, data format 4"""
    TEST_STAR = "CI-Chondrite-Lodders2019_f4.dat"
    EXPECTED_RESULT = 16373.047588932503

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
