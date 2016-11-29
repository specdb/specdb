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
    idx = igmsp.qcat.radial_search((0.038604,15.298477), 1*u.arcsec)
    assert idx >= 0
    # Blank
    idx = igmsp.qcat.radial_search((10.038604,55.298477), 1*u.arcsec)
    assert len(idx) == 0
    # Multiple (insure rank order)
    icoord = SkyCoord(ra=0.0055, dec=-1.5, unit='deg')
    idxm = igmsp.qcat.radial_search(icoord, 1*u.deg)
    coord = SkyCoord(ra=igmsp.qcat.cat['RA'][idxm], dec=igmsp.qcat.cat['DEC'][idxm], unit='deg')
    sep = icoord.separation(coord)
    isrt = np.argsort(sep)
    assert isrt[0] == 0
    assert isrt[-1] == len(idxm)-1
    # Multiple but grab only 1
    idxs = igmsp.qcat.radial_search(icoord, 1*u.deg, max=1)
    assert len(idxs) == 1


def test_match_coord(igmsp):
    # Single
    coord = SkyCoord(ra=0.038604, dec=15.298477, unit='deg')
    #
    idx = igmsp.qcat.match_coord(coord)
    assert idx[0] >= 0
    # Multiple
    coords = SkyCoord(ra=[0.038604]*2, dec=[15.298477]*2, unit='deg')
    idxs = igmsp.qcat.match_coord(coords)
    assert len(idxs) == 2

"""
def test_sdss():
    pargs = sdss_spec.parser(['751', '354'])
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    sdss_spec.main(pargs, db_file=db_file, unit_test=True)
"""
