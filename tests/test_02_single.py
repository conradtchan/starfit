import numpy as np

import starfit


def test():
    """Test a single star fits, data format 2"""
    TEST_STAR = "HE0557-4840.dat"
    EXPECTED_RESULT = 3.811527192751109

    test_result = starfit.Single(TEST_STAR, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, EXPECTED_RESULT)
