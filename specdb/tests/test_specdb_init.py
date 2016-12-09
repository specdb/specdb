# Module to run tests on instantiating specdb
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table

from specdb.specdb import SpecDB

#version = 'v01'
version = 'v02'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    # Instantiate
    sdb = SpecDB(db_file=db_file)
    assert 'BOSS_DR12' in sdb.groups

def test_load_group():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    # Instantiate
    sdb = SpecDB(db_file=db_file)
    # Load a DB
    sdb['BOSS_DR12']
    assert 'BOSS_DR12' in sdb._gdict.keys()
    # Try meta
    assert isinstance(sdb['BOSS_DR12'].meta, Table)


