import numpy as np
import pytest

import starfit

testdata = [
    ("HE1327-2326.dat", 3.292526847047224),  # data format 1
    ("HE0557-4840.dat", 3.811527192751109),  # data format 2
    ("BD_80_245.dat", 1.12390641724435),  # data format 3
    ("CI-Chondrite-Lodders2019_f4.dat", 16373.047588932503),  # data format 4
    ("CI-Chondrite-Lodders2019_f5.dat", 16372.053997257417),  # data format 5
    ("SDSS-J102915+172927.dat", 0.7040288685922813),  # data format 6
    ("CI-Chondrite-Lodders2019_f7.dat", 15616.41737463582),  # data format 7
]


@pytest.mark.parametrize("filename, expected", testdata)
def test(filename, expected):
    """Test single star fits"""

    test_result = starfit.Single(filename, silent=True).sorted_fitness[0]

    assert np.isclose(test_result, expected)
