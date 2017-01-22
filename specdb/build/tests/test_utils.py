# Module to run tests on ID methods
from __future__ import print_function, absolute_import, division, unicode_literals

import pytest
import numpy as np
import os

from astropy.coordinates import SkyCoord
from astropy.table import Table
from specdb.build import utils as spbu

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def tdata_path(filename):
    tdata_dir = os.path.join(os.path.dirname(__file__), '../../tests/files')
    return os.path.join(tdata_dir, filename)

version = 'v02'
@pytest.fixture
def igmsp():
    from specdb.specdb import IgmSpec
    db_file = tdata_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    return igmsp


def test_chk_vstack(igmsp):
    assert spbu.chk_vstack(igmsp.hdf)


def test_clean_meta():
    meta = Table()
    meta['SPEC_FILE'] = ['abas', 'adfasd', '123as']
    assert 'unicode' in meta['SPEC_FILE'].dtype.name
    # Clean
    spbu.clean_table_for_hdf(meta)
    assert 'str' in meta['SPEC_FILE'].dtype.name

def test_get_newids():
    # Faux maindb
    maindb = Table()
    maindb['RA'] = [1., 2., 3.]
    maindb['DEC'] = [2., 3, 4.]
    maindb['ID_KEY'] = np.arange(3).astype(int)
    # Input table
    meta = Table()
    meta['RA_GROUP'] = [1., 4., 6., 6., 7, 7, 2.]  # Dups but no pairs
    meta['DEC_GROUP'] = [2., 5, 7., 7., 8., 8., 3.]
    IDs = spbu.get_new_ids(maindb, meta, 'ID_KEY')
    assert IDs[0] == 0
    assert IDs[2] == IDs[2]
    assert IDs[-1] == -1
    # Fail (close pairs with a match to original)
    meta = Table()
    meta['RA_GROUP'] = [1., 4., 6., 6., 6.]
    meta['DEC_GROUP'] = [2., 5, 7., 7., 7.0002]  # Last entry adds a 'pair'
    # Fail first (close_pairs not set)
    with pytest.raises(IOError):
        IDs2 = spbu.get_new_ids(maindb, meta, 'ID_KEY')
    # Pair
    meta = Table()
    meta['RA_GROUP'] = [1., 4., 3., 3.]
    meta['DEC_GROUP'] = [2., 5, 4., 4.0002]  # Last entry adds a 'pair'
    IDs3 = spbu.get_new_ids(maindb, meta, 'ID_KEY', close_pairs=True)
    assert IDs3[-1] == 4
    # Not quite a pair (0.36" separation)
    meta = Table()
    meta['RA_GROUP'] = [1., 4., 3., 3.]
    meta['DEC_GROUP'] = [2., 5, 4., 4.0001]  # Last entry adds a 'pair'
    IDs3b = spbu.get_new_ids(maindb, meta, 'ID_KEY', close_pairs=True)
    # Now run properly
    meta = Table()
    meta['RA_GROUP'] = [1., 4., 6., 6., 6.]
    meta['DEC_GROUP'] = [2., 5, 7., 7., 7.0002]  # Last entry adds a 'pair'
    IDs = spbu.get_new_ids(maindb, meta, 'ID_KEY', close_pairs=True)





