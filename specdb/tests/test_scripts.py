# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os

from ..scripts import plot_specdb
from ..scripts import sdss_spec

#version = 'v01'
version = 'v02'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_plot_spec():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    # Various inputs
    pargs = plot_specdb.parser(['000000.45+174625.4', 'priv', '--db_file={:s}'.format(db_file)])
    plot_specdb.main(pargs, unit_test=True)
    #
    pargs = plot_specdb.parser(['00:00:00.45,+17:46:25.4', 'priv', '--db_file={:s}'.format(db_file)])
    plot_specdb.main(pargs, unit_test=True)
    #
    pargs = plot_specdb.parser(['0.00189828518376,17.7737391299', 'priv', '--db_file={:s}'.format(db_file)])
    plot_specdb.main(pargs, unit_test=True)


def test_sdss():
    pargs = sdss_spec.parser(['751', '354', 'priv'])
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    sdss_spec.main(pargs, db_file=db_file, unit_test=True)
