# Module to run tests on scripts
from __future__ import print_function, absolute_import, division, unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord

from ..specdb import IgmSpec


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def igmsp():
    version = 'v02'
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    return igmsp


def test_allspec_from_coord(igmsp):
    # One match
    spec_list, meta_list = igmsp.allspec_at_coord((0.0019, 17.7737))
    assert len(spec_list) == 1
    assert meta_list[0]['PLATE'][0] == 6173
    # Multiple matches and spectra
    spec_list, meta_list = igmsp.allspec_at_coord('001115.23+144601.8', groups=['GGG']) # Not in debug file for BOSS or SDSS
    assert spec_list[0].nspec == 2


def test_coords_to_spec(igmsp):
    coords = SkyCoord(ra=[0.0028, 0.0019], dec=[14.9747, 17.77374], unit='deg')
    spec, meta = igmsp.coords_to_spectra(coords, 'BOSS_DR12')
    # Test
    assert spec.nspec == 2
    assert meta['PLATE'][0] == 6177
    # Multiple spectra per group
    coords = SkyCoord(ra=[2.8135, 16.5802], dec=[14.7672, 0.8065], unit='deg')
    spec, meta = igmsp.coords_to_spectra(coords, 'GGG', all_spec=True)
    assert spec.nspec == 4
