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

def test_query_dict(igmsp):
    # BITWISE-OR
    qdict = {'zem': (3.,5.), 'flag_group-BITWISE-OR': [2,4,8], 'STYPE': 'QSO'}
    matches, sub_cat, IDs = igmsp.qcat.query_dict(qdict)
    assert 123 in IDs
    # BITWISE-AND
    qdict2 = {'zem': (3.,5.), 'flag_group-BITWISE-AND': [4,8], 'STYPE': 'QSO'}
    matches, sub_cat, IDs2 = igmsp.qcat.query_dict(qdict2)
    assert 123 not in IDs2
    assert 24295 in IDs2

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
    idxs = igmsp.qcat.radial_search(icoord, 1*u.deg, mt_max=1)
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


def test_cat_from_coord(igmsp):
    # Single
    coord = SkyCoord(ra=0.0019, dec=17.7737, unit='deg')
    #
    ccat = igmsp.qcat.cat_from_coords(coord)
    assert ccat['IGM_ID'][0] == 0
    # Multiple
    coords = SkyCoord(ra=[0.0028,0.0019], dec=[14.9747,17.7737], unit='deg')
    ccat2 = igmsp.qcat.cat_from_coords(coords)
    assert len(ccat2) == 2
    assert ccat2['IGM_ID'][0] == 1
    # One miss
    coords3 = SkyCoord(ra=[9.0028,0.0019], dec=[-14.9747,17.7737], unit='deg')
    ccat3 = igmsp.qcat.cat_from_coords(coords3)
    assert ccat3['IGM_ID'][0] == -1


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


