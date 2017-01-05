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



