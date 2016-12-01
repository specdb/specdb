# Module to run tests on scripts
from __future__ import print_function, absolute_import, division, unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np
import h5py

from astropy import units as u
from astropy.table import Table

from specdb.interface_group import InterfaceGroup

#version = 'v01'
version = 'v02'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_load():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    # Instantiate
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12')
    assert isinstance(boss_group.meta, Table)


