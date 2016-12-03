# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from specdb import cat_utils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_match_ids():
    tbl_IDs = np.arange(12).astype(int)
    # Simple match
    IDs = np.array([1,5,2,11])
    mIDs = cat_utils.match_ids(IDs, tbl_IDs)
    assert np.sum(mIDs < 0) == 0
    assert mIDs[0] == 1
    # Missing one when required
    IDs = np.array([1,5,2,11,99])
    with pytest.raises(IOError):
        mIDs = cat_utils.match_ids(IDs, tbl_IDs)
    # Missing when not required
    IDs = np.array([1,5,2,11,99])
    mIDs = cat_utils.match_ids(IDs, tbl_IDs, require_in_match=False)
    assert mIDs[-1] == -1
    assert mIDs[1] == 5



