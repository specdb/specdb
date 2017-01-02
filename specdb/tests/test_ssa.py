# Module to run tests on scripts
from __future__ import print_function, absolute_import, division, unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os
import numpy as np

from astropy.io.votable.tree import VOTableFile

from ..specdb import IgmSpec
from specdb import ssa as spdb_ssa


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def igmsp():
    version = 'v02'
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    igmsp = IgmSpec(db_file=db_file)
    return igmsp


def test_empty_vo():
    evot = spdb_ssa.empty_vo()
    assert evot.resources[0].type == 'results'
    assert len(evot.resources[0].infos) == 0
    assert len(evot.resources[0].params) == 0
    #
    evot = spdb_ssa.empty_vo(rtype='meta')
    assert evot.resources[0].type == 'meta'


def test_input_params():
    evot = spdb_ssa.empty_vo()
    iparams = spdb_ssa.input_params(evot)
    names = [param.name for param in iparams]
    # Test
    assert 'INPUT:POS' in names
    assert 'INPUT:SIZE' in names
    assert 'INPUT:BAND' in names
    assert 'INPUT:FORMAT' in names


def test_metaquery_params():
    oparams, pIDs = spdb_ssa.metaquery_param()
    #
    assert 'SpatialLocation' in pIDs


def test_default_fields():
    title = 'BOSS: DR12 Quasars'
    ssa_dict = spdb_ssa.default_fields(title)
    assert ssa_dict['FluxUcd'] == 'arith.ratio;phot.flux.density'
    assert ssa_dict['FluxUnit'] == ''
    assert ssa_dict['SpecUcd'] == 'em.wl'
    assert ssa_dict['SpecUnit'] == 'Angstrom'
    # flambda
    ssa_dict = spdb_ssa.default_fields(title, flux='flambda')
    assert ssa_dict['FluxUcd'] == 'phot.fluDens;em.wl'
    # calib
    ssa_dict = spdb_ssa.default_fields(title, flux='flambda', fxcalib='ABSOLUTE')
    assert ssa_dict['FluxCalib'] == 'ABSOLUTE'


def test_ssa_init(igmsp):
    ssai = spdb_ssa.SSAInterface(igmsp)
    assert ssai.name == 'SSAI_igmspec'


def test_query_data(igmsp):
    ssai = spdb_ssa.SSAInterface(igmsp)
    votable = ssai.querydata('0.0019,17.7737', SIZE=1e-3)
    # Test
    assert isinstance(votable, VOTableFile)
    assert len(votable.resources[0].tables) == 1
    votable.to_xml(data_path('tst.xml'))


def test_metadata(igmsp):
    ssai = spdb_ssa.SSAInterface(igmsp)
    votable = ssai.querydata(FORMAT='METADATA')
    # Test
    assert isinstance(votable, VOTableFile)
    assert len(votable.resources[0].tables) == 0

