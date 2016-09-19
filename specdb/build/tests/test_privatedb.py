# Module to run tests on scripts

import pytest
import numpy as np
import os
import h5py

from specdb.privatedb import build as pbuild


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_grab_files():
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles = pbuild.grab_files(data_dir)
    #
    assert len(ffiles) == 2


def test_meta():
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles = pbuild.grab_files(data_dir)
    meta = pbuild.mk_meta(ffiles, fname=True, skip_badz=True)
    #
    np.testing.assert_allclose(meta['zem'].data, (2.39499998093, 2.59719920158))


def test_ingest():
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles = pbuild.grab_files(data_dir)
    meta = pbuild.mk_meta(ffiles, fname=True, skip_badz=True)
    hdf = h5py.File('tmp.hdf5','w')
    pbuild.ingest_spectra(hdf, 'test', meta)
    hdf.close()
    # Read
    tmp = h5py.File('tmp.hdf5','r')
    # Test
    assert 'meta' in tmp['test'].keys()
    assert isinstance(tmp['test/spec'].value, np.ndarray)

