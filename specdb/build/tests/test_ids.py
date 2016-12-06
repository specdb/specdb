# Module to run tests on ID methods

import pytest
import numpy as np
import os

from astropy.coordinates import SkyCoord
from astropy.table import Table
from specdb.build import utils as spbu


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_get_newids():
    # Faux maindb
    maindb = Table()
    maindb['RA'] = [1., 2., 3.]
    maindb['DEC'] = [2., 3, 4.]
    maindb['ID_KEY'] = np.arange(3).astype(int)
    # Input table
    meta = Table()
    meta['RA_SPEC'] = [1., 4., 6., 6., 7, 7, 2.]  # Dups but no pairs
    meta['DEC_SPEC'] = [2., 5, 7., 7., 8., 8., 3.]
    IDs = spbu.get_new_ids(maindb, meta, 'ID_KEY')
    assert IDs[0] == 0
    assert IDs[2] == IDs[2]
    assert IDs[-1] == -1

    # Pair
    meta = Table()
    meta['RA_SPEC'] = [1., 4., 6., 6., 6.]
    meta['DEC_SPEC'] = [2., 5, 7., 7., 7.0001]  # Last entry adds a 'pair'
    # Fail first (close_pairs not set)
    with pytest.raises(IOError):
        IDs2 = spbu.get_new_ids(maindb, meta, 'ID_KEY')
    # Fail (close pairs with a match to original)
    meta = Table()
    meta['RA_SPEC'] = [1., 4., 3., 3.]
    meta['DEC_SPEC'] = [2., 5, 4., 4.0001]  # Last entry adds a 'pair'
    IDs3 = spbu.get_new_ids(maindb, meta, 'ID_KEY', close_pairs=True)
    # Now run properly
    meta = Table()
    meta['RA_SPEC'] = [1., 4., 6., 6., 6.]
    meta['DEC_SPEC'] = [2., 5, 7., 7., 7.0002]  # Last entry adds a 'pair'
    IDs = spbu.get_new_ids(maindb, meta, 'ID_KEY', close_pairs=True)
    pytest.set_trace()





