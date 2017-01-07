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


def test_init():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    # Instantiate
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12', 'IGM_ID')
    assert isinstance(boss_group.meta, Table)
    hdf.close()

def test_getrows_one_source():
    # One source for 1 row
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12', 'IGM_ID')
    ID = 0
    rows = boss_group.ids_to_allrows(ID)
    assert rows.size == 1
    # One source for multiple spectra
    ggg_group = InterfaceGroup(hdf, 'GGG', 'IGM_ID')
    ID = 17656
    rows = ggg_group.ids_to_allrows(ID)
    assert rows.size == 2
    # One source, first spectrum
    rows = ggg_group.ids_to_firstrow(ID)
    assert rows.size == 1
    hdf.close()


def test_getrows_sources():
    # Init
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12', 'IGM_ID')
    # All rows, unique
    IDs = np.array([0,2])
    rows = boss_group.ids_to_allrows(IDs)
    assert rows.size == 2
    # Dealing with multiple spectra
    ggg_group = InterfaceGroup(hdf, 'GGG', 'IGM_ID')
    IDs = np.array([32720,17656])
    rows = ggg_group.ids_to_allrows(IDs)
    assert rows.size == 4
    rows = ggg_group.ids_to_firstrow(IDs)
    assert rows.size == 2
    assert rows[0] > rows[1]
    hdf.close()


def test_grab_specmeta():  # Use rows only
    # Init
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12', 'IGM_ID')
    # Get spec and meta from one row
    rows = 1
    spec, meta = boss_group.grab_specmeta(rows)
    assert spec.nspec == 1
    # Get spec and meta from multiple rows, with repeat
    rows = np.array([4,3,3,1,4])
    spec, meta = boss_group.grab_specmeta(rows)
    assert spec.nspec == 5
    assert meta['IGM_ID'][0] == meta['IGM_ID'][-1]


def test_spec_from_meta():
    # Init
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    hdf = h5py.File(db_file, 'r')
    boss_group = InterfaceGroup(hdf, 'BOSS_DR12', 'IGM_ID')
    meta = boss_group.meta
    # One entry
    imeta = meta[0:1]
    spec = boss_group.spec_from_meta(imeta)
    assert spec.nspec == 1
    # Get spec and meta from multiple rows, with repeat
    rows = np.array([4,3,3,1,4])
    imeta = meta[rows]
    spec = boss_group.spec_from_meta(imeta)
    assert spec.nspec == 5
