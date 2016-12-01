""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import h5py
import pdb

from astropy.table import Table

from specdb import defs


def match_ids(IDs, tbl_IDs):
    """ Match input IDs to another array of IDs (usually in a table)

    Parameters
    ----------
    IDs : ndarray
    tbl_IDs : ndarray

    Returns
    -------

    """
    rows = -1 * np.ones_like(IDs).astype(int)
    # Find which IDs are in tbl_IDs
    in_tbl = np.in1d(IDs, tbl_IDs)
    rows[~in_tbl] = -1
    #
    IDs_intbl = IDs[in_tbl]
    # Find indices of input IDs in meta table -- first instance in meta only!
    xsorted = np.argsort(tbl_IDs)
    ypos = np.searchsorted(tbl_IDs, IDs_intbl, sorter=xsorted)
    indices = xsorted[ypos]
    rows[in_tbl] = indices
    return rows

def flag_to_surveys(flag, survey_dict):
    """ Convert flag_survey to list of surveys

    Parameters
    ----------
    flag : int
    survey_dict : dict

    Returns
    -------
    surveys : list

    """
    surveys = []
    for key,sflag in survey_dict.items():
        if flag % (2*sflag) >= sflag:
            surveys.append(key)
    # Return
    return surveys


def write_cat_to_fits(DB_file, cat_fits_file):
    """ Simple script to write the catalog file to a FITS file (mainly for others)
    Parameters
    ----------
    DB_file : str
      Full path to the DB file which contains the catalog
    cat_fits_file : str
      Filename for the FITS file

    Returns
    -------

    """
    if '.fits' not in cat_fits_file:
        raise IOError("Output file {:s} must have .fits extension".format(cat_fits_file))
    # Read
    hdf = h5py.File(DB_file, 'r')
    cat = Table(hdf['catalog'])
    # Write
    cat.write(cat_fits_file)
    # Finish
    hdf.close()
    return
