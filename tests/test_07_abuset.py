import glob

import pytest

from starfit.autils.abuset import AbuSet

testdata = [data for data in glob.glob("src/starfit/data/ref/*.dat")]


@pytest.mark.parametrize("filename", testdata)
def test(filename):
    """Test loading of all star data"""

    AbuSet(filename)
