# Module to run tests on meta query
from __future__ import print_function, absolute_import, division, unicode_literals

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
    # GGG
    qdict = {'TELESCOPE': 'Gemini-North', 'NPIX': (1580,1583), 'GRATING': ['B600', 'R400']}
    qmeta = igmsp.query_meta(qdict)
    assert qmeta['GROUP'][0] == 'GGG'


def test_meta_from_position(igmsp):
    # One match
    meta = igmsp.meta_from_position((0.0019,17.7737), 1*u.arcsec)
    assert len(meta) == 1
    # Blank
    meta2 = igmsp.meta_from_position((10.038604,55.298477), 1*u.arcsec)
    assert meta2 is None
    # Multiple sources (insure rank order)
    meta3 = igmsp.meta_from_position((0.0055,-1.5), 1*u.deg)
    assert len(meta3) == 2
    assert np.isclose(meta3['R'][0],meta3['R'][1])
    # Multiple meta entries (GGG)
    meta4 = igmsp.meta_from_position('001115.23+144601.8', 1*u.arcsec)
    assert len(meta4) == 2
    assert meta4['R'][0] != meta4['R'][1]
    # Multiple but grab closest source
    meta5 = igmsp.meta_from_position((0.0055,-1.5), 1*u.deg, max_match=1)
    assert len(meta5) == 1
    # Groups
    meta = igmsp.meta_from_position((2.813500,14.767200), 20*u.deg, groups=['GGG','HD-LLS_DR1'])
    for group in meta['GROUP'].data:
        assert group in ['GGG', 'HD-LLS_DR1']


def test_meta_from_coords(igmsp):
    # Single
    coord = SkyCoord(ra=0.0019, dec=17.7737, unit='deg')
    _, meta = igmsp.meta_from_coords(coord)
    meta['GROUP'] == 'BOSS_DR12'
    # Single miss
    coord = SkyCoord(ra=0.0019, dec=-17.7737, unit='deg')
    matches, meta = igmsp.meta_from_coords(coord)
    assert meta is None
    assert np.sum(matches) == 0
    # Multiple sources
    coords = SkyCoord(ra=[0.0028,0.0019], dec=[14.9747,17.7737], unit='deg')
    _, meta = igmsp.meta_from_coords(coords)
    assert len(meta) == 2
    # With one query retrieving None
    matchesN, metaN = igmsp.meta_from_coords(coords, query_dict=dict(PLATE=6177))
    assert np.sum(matchesN) == 1
    assert metaN['RA_GROUP'][1].mask == True
    # Multiple sources with one bad
    coords = SkyCoord(ra=[0.0028,0.0019], dec=[-14.9747,17.7737], unit='deg')
    matches, meta = igmsp.meta_from_coords(coords)
    assert matches[0] == False
    # Multiple hits with first=True
    coord = SkyCoord(ra=2.813458, dec=14.767167, unit='deg')
    _, meta4 = igmsp.meta_from_coords(coord)
    assert len(meta4) == 1  # Only grabs first
    assert meta4['GROUP_ID'] == 0
    #
    coords = SkyCoord(ra=[0.0028,9.99,2.813458], dec=[14.9747,-9.99,14.767167], unit='deg')
    matches, meta4b = igmsp.meta_from_coords(coords)
    assert meta4b['IGM_ID'][1] == -1
    assert meta4b['RA_GROUP'][1].mask == True
    #
    coords = SkyCoord(ra=[2.813458]*2, dec=[14.767167]*2, unit='deg')
    _, meta5 = igmsp.meta_from_coords(coords)
    assert len(meta5) == 2
    assert np.sum(meta5['GROUP_ID'] == meta5['GROUP_ID']) == 2
    # Multiple hits on single source with first=False
    coord = SkyCoord(ra=2.813458, dec=14.767167, unit='deg')
    _, meta6 = igmsp.meta_from_coords(coord, first=False)
    assert len(meta6) == 1
    assert len(meta6[0]) == 2
    assert meta6[0]['GROUP_ID'][0] == 0
    # Multiple hits on two sources with first=False
    coords = SkyCoord(ra=[0.0028,2.813458], dec=[14.9747,14.767167], unit='deg')
    matches7, meta7 = igmsp.meta_from_coords(coords, first=False)
    assert len(meta7[0]) == 1
    assert len(meta7[1]) == 2
    assert np.sum(matches7) == 2
    # Multiple hits on two sources with first=False; limit by groups
    coords = SkyCoord(ra=[0.0028,2.813458], dec=[14.9747,14.767167], unit='deg')
    matches7b, meta7b = igmsp.meta_from_coords(coords, first=False, groups=['GGG'])
    assert meta7b[0] is None
    # Multiple hits on mixed sources with first=False
    coords = SkyCoord(ra=[0.0028,9.99,2.813458], dec=[14.9747,-9.99,14.767167], unit='deg')
    matches8, meta8 = igmsp.meta_from_coords(coords, first=False)
    assert meta8[1] is None
    assert np.sum(matches8) == 2
    # Limit by groups
    matches8b, meta8b = igmsp.meta_from_coords(coords, first=False, groups=['GGG'])
