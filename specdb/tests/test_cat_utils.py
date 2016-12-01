# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np
import h5py

from astropy import units as u
from astropy.table import Table

from specdb import cat_utils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_match_ids():
    IDs = np.array([1,5,2,11])
    tbl_IDs = 3 + np.arange(10).astype(int)
    # Match
    mIDs = cat_utils.match_ids(IDs, tbl_IDs)
    assert mIDs[0] == -1
    assert mIDs[1] == 2



