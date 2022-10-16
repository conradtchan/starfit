import numpy as np

import starfit


def test():
    """Test a single star fits, data format 6"""
    TEST_STAR = "SDSS-J102915+172927.dat"
    EXPECTED_RESULT = 0.7040288685922813

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
