# Module to run tests on scripts

import pytest
import numpy as np
import os
import h5py

from astropy.table import Table
from specdb.build import privatedb as pbuild


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_grab_files():
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles, _ = pbuild.grab_files(data_dir)
    #
    assert len(ffiles) == 2


def test_meta():
    ztbl = Table.read(os.path.join(os.path.dirname(__file__), 'files', 'ztbl_E.fits'))
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles,_ = pbuild.grab_files(data_dir)
    meta = pbuild.mk_meta(ffiles, ztbl, fname=True, skip_badz=True, mdict=dict(INSTR='HIRES'))
    #
    np.testing.assert_allclose(meta['zem'].data, (2.39499998093, 2.59719920158))


def test_ingest():
    ztbl = Table.read(os.path.join(os.path.dirname(__file__), 'files', 'ztbl_E.fits'))
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    ffiles,_ = pbuild.grab_files(data_dir)
    meta = pbuild.mk_meta(ffiles, ztbl, fname=True, skip_badz=True, mdict=dict(INSTR='HIRES'))
    hdf = h5py.File('tmp.hdf5','w')
    pbuild.ingest_spectra(hdf, 'test', meta)
    hdf.close()
    # Read
    tmp = h5py.File('tmp.hdf5','r')
    # Test
    assert 'meta' in tmp['test'].keys()
    assert isinstance(tmp['test/spec'].value, np.ndarray)


def test_mkdb():
    import specdb
    from linetools import utils as ltu
    # Generate a fake ztbl
    obj = ['J220758.30+125944.3',
           'J172524.66+303803.9',
           'J230044.36+015541.7',
           'J095240.17+515250.03',
           'J095243.05+515121.15',
           ]
    ztbl = Table()
    RA, DEC, ZEM, ZSOURCE = [], [], [], []
    for kk,iobj in enumerate(obj):
           icoord = ltu.radec_to_coord(iobj)
           RA.append(icoord.ra.value)
           DEC.append(icoord.dec.value)
           ZEM.append(1 + kk*0.1)
           ZSOURCE.append('UNKNW')
    ztbl['RA'] = RA
    ztbl['DEC'] = DEC
    ztbl['ZEM'] = ZEM
    ztbl['ZEM_SOURCE'] = ZSOURCE

    # Run
    tree = specdb.__path__[0]+'/data/privateDB'
    pbuild.mk_db('tst_db', tree, 'tst_db.hdf5', ztbl, fname=True)
