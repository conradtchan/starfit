import glob

import pytest

from starfit.autils.stardb import StarDB

testdata = [data for data in glob.glob("src/starfit/data/db/*.stardb.*")]


@pytest.mark.parametrize("filename", testdata)
def test(filename):
    """Test loading of all star data"""

    StarDB(filename)
