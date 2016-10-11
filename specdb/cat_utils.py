""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import h5py
import pdb

from astropy.table import Table

from specdb import defs

def flag_to_surveys(flag, survey_dict=None):
    """ Convert flag_survey to list of surveys

    Parameters
    ----------
    flag : int

    Returns
    -------
    surveys : list

    """
    if survey_dict is None:
        survey_dict = defs.get_survey_dict()
    #
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
