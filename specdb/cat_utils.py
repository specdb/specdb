""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import h5py
import pdb

from astropy.table import Table


def match_ids(IDs, match_IDs, require_in_match=True):
    """ Match input IDs to another array of IDs (usually in a table)
    Return the rows aligned with input IDs

    Parameters
    ----------
    IDs : ndarray
    match_IDs : ndarray
    require_in_match : bool, optional
      Require that each of the input IDs occurs within the match_IDs

    Returns
    -------
    rows : ndarray
      Rows in match_IDs that match to IDs, aligned
      -1 if there is no match

    """
    rows = -1 * np.ones_like(IDs).astype(int)
    # Find which IDs are in match_IDs
    in_match = np.in1d(IDs, match_IDs)
    if require_in_match:
        if np.sum(~in_match) > 0:
            raise IOError("qcat.match_ids: One or more input IDs not in match_IDs")
    rows[~in_match] = -1
    #
    IDs_inmatch = IDs[in_match]
    # Find indices of input IDs in meta table -- first instance in meta only!
    xsorted = np.argsort(match_IDs)
    ypos = np.searchsorted(match_IDs, IDs_inmatch, sorter=xsorted)
    indices = xsorted[ypos]
    rows[in_match] = indices
    return rows


def flags_to_groups(flags, group_dict):
    """ Convert flag values to np.array of group names comma separated

    Parameters
    ----------
    flag : list or ndarray
      BITWISE flags
    group_dict : dict
      dict that converts the BITWISE flags to group names

    Returns
    -------
    groups : str ndarray
      Array of expanded list of groups corresponding to input flags
    """
    # Ugly for loops appear necessary
    group_list = [[] for ii in range(len(flags))]
    for key,value in group_dict.items():
        in_flags = np.where(flags & value)[0]
        for ii in in_flags:
            group_list[ii].append(key)
    # More
    for ii,ilist in enumerate(group_list):
        group_list[ii] = ','.join(ilist)
    # Finally
    garray = np.array(group_list)
    # Return
    return garray


def flag_to_groups(flag, group_dict):
    """ Convert flag value to list of groups

    Parameters
    ----------
    flag : int
      BITWISE flag
    group_dict : dict
      dict that converts the BITWISE flags

    Returns
    -------
    groups : list
      Expanded list of groups corresponding to input flag
    """
    groups = []
    for key,sflag in groups.items():
        if flag % (2*sflag) >= sflag:
            groups.append(key)
    # Return
    return groups


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
