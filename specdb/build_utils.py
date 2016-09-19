""" Module to build the hdf5 database file for ExGalSpec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import igmspec
import os

import h5py
import numbers, json
import pdb

from igmspec import defs
from igmspec.ingest import boss, hdlls, kodiaq, ggg, sdss, hst_z2, myers, twodf, xq100
from igmspec.ingest import hdla100
from igmspec.ingest import esidla
from igmspec.ingest import cos_halos

from astropy.table import Table, vstack, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

from linetools import utils as ltu


