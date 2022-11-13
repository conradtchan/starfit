import glob

import pytest

from starfit.star import Star

testdata = [data for data in glob.glob("src/starfit/data/stars/*.dat")]


@pytest.mark.parametrize("filename", testdata)
def test(filename):
    """Test loading of all star data"""

    Star(filename)
