"""
Definition of data and files for the tests.
The files are located in the data directory.
"""
import os
from os.path import join as pjoin

test_dir = os.path.dirname(os.path.abspath(__file__))  # directory of test files
data_dir = pjoin(test_dir, 'data')  # directory of data for tests


MICROARRAY_FLUTYPE_TEST_DIR = pjoin(data_dir, 'microarray/flutype_test')
MICROARRAY_FLUTYPE_TEST_DATA_ID = 'flutype_test'

