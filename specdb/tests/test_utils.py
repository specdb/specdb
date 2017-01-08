# Module to run tests on cat_utils
from __future__ import print_function, absolute_import, division, unicode_literals

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from specdb import utils

from ..specdb import IgmSpec

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

version = 'v02'
@pytest.fixture
def igmsp():
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    return igmsp


def test_clean_vstack(igmsp):
    def run_tst(sdb):
        # Load
        tbls = []
        for group in sdb.groups:
            tbls.append(sdb[group].meta[0:1])
        # Stack
        stack = utils.clean_vstack(tbls, sdb.groups)
    # Internal
    run_tst(igmsp)
    # Full IGMSpec
    if os.getenv('IGMSPEC_DB') is not None:
        tigmsp = IgmSpec()
        run_tst(tigmsp)

