"""
Utility and helper functions.
"""
from __future__ import print_function, absolute_import
import os.path
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
import cv2
import flutype_analysis

PATTERN_PEP_GAL = "{}_pep.gal"
PATTERN_VIR_GAL = "{}_vir.gal"
PATTERN_DATA = "{}.csv"
PATTERN_STD = "{}_std.csv"
PATTERN_TIF = "{}_600_100_635.tif"


def load_data(data_id, directory='data/flutype_test', what="all"):
    """ Loads datasets for given data_id.
    See naming schema above.

    :param data_id: id of the Data
    :param
    :return: dictionary of information
    """
    print("-" * 80)
    print("Loading data corresponding to data_id: <{}> in dir <{}>".format(data_id, directory))
    print("-" * 80)

    d = {'data_id': data_id}

    f_gal_vir = os.path.join(directory, PATTERN_VIR_GAL.format(data_id))
    d['gal_vir'] = pd.read_csv(f_gal_vir, sep='\t', index_col="ID")

    print("{}:{}".format("Virus .gal", f_gal_vir))

    f_gal_pep = os.path.join(directory, PATTERN_PEP_GAL.format(data_id))
    d['gal_pep'] = pd.read_csv(f_gal_pep, sep='\t', index_col="ID")

    print("Peptide .gal :{}".format(f_gal_pep))

    f_data = os.path.join(directory, PATTERN_DATA.format(data_id))
    if os.path.exists(f_data):
        d['data'] = pd.read_csv(f_data, sep='\t', index_col=0)
        print("Spot intensity file  :{}".format(f_data))

    f_tif = os.path.join(directory, PATTERN_TIF.format(data_id))
    if os.path.exists(f_tif):
        d['tif'] = cv2.imread(f_tif, 0)
        print("Image file  :{}".format(f_tif))

    f_std = os.path.join(directory, PATTERN_STD.format(data_id))

    if os.path.exists(f_std):
        d['data_std'] = pd.read_csv(f_std, sep='\t', index_col=0)
        print("Intensity standard deviation:{}".format(f_std))
    else:
        print("Spot intensities for the data ID ({}) are not averaged but primary values".format(data_id))

    return d

def map_strings_to_number(strings):
    """Transforms list of strings into numbers."""
    counter = 1
    d = {}

    for s in strings:
        if s in d:
            pass
        else:
            d[s] = counter
            counter += 1

    return d


def discrete_cmap(N, base_cmap=None):
    """ Create an N-bin discrete colormap from the specified input map"""
    # Note that if base_cmap is a String or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for String, None, or a colormap instance:
    if N == 1:
        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N + 1))
        cmap_name = base.name + str(N)

    else:

        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N))
        cmap_name = base.name + str(N)

    return base.from_list(cmap_name, color_list, N)

def assure_path_exists(path):
    """
    creates folder if not existing

    :param path:
    :return:
    """
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)

def merge_dataset(data1_id,data1_dir,data2_id, data2_dir, name):
    assure_path_exists("data/" + name + "/")
    data1 = load_data(data_id=data1_id, directory=data1_dir, what="all")
    data2 = load_data(data_id=data2_id, directory=data2_dir, what="all")
    ana1 = flutype_analysis.analysis.Analysis(data1)
    ana2 = flutype_analysis.analysis.Analysis(data2)

    ana1.spots()


















