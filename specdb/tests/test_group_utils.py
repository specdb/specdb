# Module to run tests on group_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.table import Table

from specdb import group_utils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_show_meta():
    # Dummy meta
    meta = Table()
    meta['IGM_ID'] = [134,2234,211]
    meta['GROUP_ID'] = [0,1,2]
    meta['RA'] = [122.2234, 200.22222, 9.938444]
    meta['DEC'] = [-13.22333, 23.211199, 55.2323232]
    meta['zem'] = [1.23244, 2.220848, 10.188888]
    meta['SPEC_FILE'] = ['J1.fits', 'J2.fits', 'J3.fits']
    # Show
    group_utils.show_group_meta(meta)
    # Now with keys
    group_utils.show_group_meta(meta, meta_keys=['IGM_ID', 'RA', 'zem'])

