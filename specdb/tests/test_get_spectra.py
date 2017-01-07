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


def test_spectra_from_meta(igmsp):  # Base level operation
    # One match
    meta = igmsp.meta_from_position((0.0019,17.7737), 1*u.arcsec)
    spec = igmsp.spectra_from_meta(meta)
    # Two sources, one meta entry each
    meta2 = igmsp.meta_from_position((0.0055,-1.5), 1*u.deg)
    spec2 = igmsp.spectra_from_meta(meta2)
    assert spec2.nspec == 2
    # Many sources and meta entries;  groups separated
    meta3 = igmsp.meta_from_position((2.813500,14.767200), 20*u.deg)#, groups=['GGG','HD-LLS_DR1'])
    spec3 = igmsp.spectra_from_meta(meta3)
    assert spec3.nspec == 15
    # Many sources and meta entries;  groups scrambled
    idx = np.arange(15).astype(int)
    idx[1] = 13
    idx[13] = 1
    meta4 = meta3[idx]
    spec4 = igmsp.spectra_from_meta(meta4)#, debug=True)
    spec4.select = 1
    assert np.isclose(meta4['WV_MIN'][1], spec4.wvmin.value)


def test_spectra_from_coord(igmsp):
    # No match
    specN, metaN = igmsp.spectra_from_coord((0.0019, -17.7737))
    assert specN is None
    assert metaN is None
    # One match
    spec, meta = igmsp.spectra_from_coord((0.0019, 17.7737))
    assert spec.nspec == 1
    assert meta['PLATE'][0] == 6173
    # Multiple matches and spectra
    spec, meta = igmsp.spectra_from_coord('001115.23+144601.8')#, groups=['GGG']) # Not in debug file for BOSS or SDSS
    assert spec.nspec == 2
    # Many matches; takes closest
    spec, meta = igmsp.spectra_from_coord((0.0019, 17.7737), tol=10*u.deg)
    assert spec.nspec == 1
    assert meta['PLATE'][0] == 6173


def test_spectra_from_ID(igmsp):
    spec, meta = igmsp.spectra_from_ID(3244)
    assert spec.nspec == 2


def test_spectra_in_group(igmsp):
    # Missed a source -- raises IOError
    coords = SkyCoord(ra=[0.0028, 0.0019], dec=[14.9747, -17.77374], unit='deg')
    with pytest.raises(IOError):
        spec, meta = igmsp.spectra_in_group(coords, 'BOSS_DR12')
    # Another with both missing
    coords = SkyCoord(ra=[2.8135,16.5802], dec=[-14.7672, -0.8065], unit='deg')
    with pytest.raises(IOError):
        spec, meta = igmsp.spectra_in_group(coords, 'GGG')
    # Each source has only spectrum in the group
    coords = SkyCoord(ra=[0.0028, 0.0019], dec=[14.9747, 17.77374], unit='deg')
    spec, meta = igmsp.spectra_in_group(coords, 'BOSS_DR12')
    # Test
    assert spec.nspec == 2
    assert meta['PLATE'][0] == 6177
    # Each source has multiple spectra in the group
    coords = SkyCoord(ra=[2.8135,16.5802], dec=[14.7672, 0.8065], unit='deg')
    spec, meta = igmsp.spectra_in_group(coords, 'GGG')
    assert meta['DISPERSER'][0] == 'B600'
    qdict = dict(DISPERSER='R400')
    spec, meta = igmsp.spectra_in_group(coords, 'GGG', query_dict=qdict)
    assert meta['DISPERSER'][0] == 'R400'
    # Another with bad grating
    qdict = dict(DISPERSER='X400')
    with pytest.raises(IOError):
        spec, meta = igmsp.spectra_in_group(coords, 'GGG', query_dict=qdict)
    '''
    # Multiple spectra per group
    coords = SkyCoord(ra=[2.8135, 16.5802], dec=[14.7672, 0.8065], unit='deg')
    spec, meta = igmsp.coords_to_spectra(coords, 'GGG', all_spec=True)
    assert spec.nspec == 4
    '''
