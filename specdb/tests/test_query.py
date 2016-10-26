# Module to run tests on scripts
from __future__ import print_function, absolute_import, division, unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os

from astropy import units as u
from astropy.coordinates import SkyCoord

from ..specdb import IgmSpec

#version = 'v01'
version = 'v02'


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_radial_search():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    idx = igmsp.qcat.radial_search((0.038604,15.298477), 1*u.arcsec)
    assert idx >= 0
    # Blank
    idx = igmsp.qcat.radial_search((10.038604,55.298477), 1*u.arcsec)
    assert len(idx) == 0


def test_match_coord():
    # Single
    coord = SkyCoord(ra=0.038604, dec=15.298477, unit='deg')
    #
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
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
