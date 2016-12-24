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

#version = 'v01'
version = 'v02'


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def igmsp():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    return igmsp

def test_radial_search(igmsp):
    # One match
    idx = igmsp.qcat.radial_search((0.0019,17.7737), 1*u.arcsec)
    assert idx >= 0
    # Blank
    idx = igmsp.qcat.radial_search((10.038604,55.298477), 1*u.arcsec)
    assert len(idx) == 0
    # Multiple (insure rank order)
    icoord = SkyCoord(ra=0.0055, dec=-1.5, unit='deg')
    idxm = igmsp.qcat.radial_search(icoord, 1*u.deg)
    # Test
    ras = []
    decs = []
    for ii in idxm:
        mt = np.where(igmsp.cat['IGM_ID'] == ii)[0]
        ras.append(igmsp.cat['RA'][mt][0])
        decs.append(igmsp.cat['DEC'][mt][0])
    coord = SkyCoord(ra=ras, dec=decs, unit='deg')
    sep = icoord.separation(coord)
    isrt = np.argsort(sep)
    assert isrt[0] == 0
    assert isrt[-1] == len(idxm)-1
    # Multiple but grab only 1
    idxs = igmsp.qcat.radial_search(icoord, 1*u.deg, max=1)
    assert len(idxs) == 1


def test_match_coord(igmsp):
    # Single
    coord = SkyCoord(ra=0.0019, dec=17.7737, unit='deg')
    #
    idx = igmsp.qcat.match_coord(coord)
    assert idx[0] >= 0
    # Multiple
    coords = SkyCoord(ra=[0.0019]*2, dec=[17.7737]*2, unit='deg')
    idxs = igmsp.qcat.match_coord(coords)
    assert len(idxs) == 2
    # Dataset
    idxs2 = igmsp.qcat.match_coord(coords, group='BOSS_DR12')
    assert np.sum(idxs2 >= 0) == 2
    idxs3 = igmsp.qcat.match_coord(coords, group='HD-LLS_DR1')
    assert np.sum(idxs3 >= 0) == 0

def test_chk_in_group(igmsp):
    # BOSS
    answer, query = igmsp.qcat.chk_in_group(np.array([0,1,2]), 'BOSS_DR12')
    assert answer
    assert query.size == 3
    #
    answer, query = igmsp.qcat.chk_in_group(np.array([0,1,2]), 'GGG')
    assert not answer

def test_ids_in_groups(igmsp):
    # BOSS
    IDs, mask = igmsp.qcat.find_ids_in_groups(['BOSS_DR12'])
    assert IDs.size == 19
    assert IDs[0] == 0
    # BOSS and HD-LLS -- Both
    IDs2, _ = igmsp.qcat.find_ids_in_groups(['HD-LLS_DR1', 'GGG'])
    assert IDs2.size == 1
    # BOSS and HD-LLS -- Either
    IDs3, _ = igmsp.qcat.find_ids_in_groups(['HD-LLS_DR1', 'GGG'], in_all=False)
    assert IDs3.size == 9
    # With input IDs
    IDs4, _ = igmsp.qcat.find_ids_in_groups(['BOSS_DR12'], np.array([0,1]))
    assert IDs4.size == 2


