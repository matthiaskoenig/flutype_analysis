from __future__ import print_function, absolute_import

from flutype_analysis import analysis ,utils
from flutype_analysis.tests import testdata


def test_load_data():
    directory = testdata.MICROARRAY_FLUTYPE_TEST_DIR
    data_id = testdata.MICROARRAY_FLUTYPE_TEST_DATA_ID

    d = utils.load_data(data_id, directory)
    assert len(d) == 4
    assert 'data_id' in d
    assert 'data' in d
    assert 'gal_pep' in d
    assert 'gal_vir' in d


def test_create_spot():
    directory = testdata.MICROARRAY_FLUTYPE_TEST_DIR
    data_id = testdata.MICROARRAY_FLUTYPE_TEST_DATA_ID
    d = utils.load_data(data_id, directory)
    # creates the spot DataFrame
    ana = analysis.Analysis(d)

    assert hasattr(ana, 'spot')
    assert ana.spot is not None
    assert len(ana.spot) > 0


# FIXME: This is not working yet, setup and teardown for pytest needed
'''
class TestAnalysis:
    def setup_method(self, test_method):
        directory = testdata.MICROARRAY_FLUTYPE_TEST_DIR
        data_id = testdata.MICROARRAY_FLUTYPE_TEST_DATA_ID
        d = analysis.load_data(data_id, directory)
        # creates the spot DataFrame
        self.ana = analysis.Analysis(d)

    def teardown_method(self, test_method):
        self.ana = None

    def test_heatmap(self):
        fig = ana.heatmap(figsize=(20, 10))
        assert fig is not None
'''