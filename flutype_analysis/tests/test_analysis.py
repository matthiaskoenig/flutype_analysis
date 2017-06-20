from __future__ import print_function, absolute_import

from flutype_analysis import analysis
from flutype_analysis.tests import testdata


def test_load_data():
    directory = testdata.MICROARRAY_FLUTYPE_TEST_DIR
    data_id = testdata.MICROARRAY_FLUTYPE_TEST_DATA_ID

    d = analysis.load_data(data_id, directory)
    assert len(d) == 5
    assert 'data_id' in d
    assert 'data' in d
    assert 'meta' in d
    assert 'gal_pep' in d
    assert 'gal_vir' in d
