""" Module of utilities for building DB files
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

import numbers
import pdb

from igmspec import defs

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

#from linetools import utils as ltu



def add_to_flag(cur_flag, add_flag):
    """ Add a bitwise flag to an existing flat

    Parameters
    ----------
    cur_flag : int or ndarray
    add_flag : int

    Returns
    -------
    new_flag : int or ndarray

    """
    if isinstance(cur_flag, numbers.Integral):
        if (cur_flag % add_flag) < (add_flag//2):
            cur_flag += add_flag
        return cur_flag
    else:  # Array
        mods = cur_flag % add_flag
        upd = mods < (add_flag//2)
        cur_flag[upd] += add_flag
        return cur_flag


def chk_maindb_join(maindb, newdb):
    """Check that new data is consistent with existing table

    Parameters
    ----------
    maindb : Table
    newdb : Table

    Returns
    -------
    answer : bool

    """
    # One way
    for key in newdb.keys():
        try:
            assert key in maindb.keys()
        except AssertionError:
            pdb.set_trace()
            return False
    return True


def chk_for_duplicates(maindb):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table

    Return
    ------
    result : bool
      * True = pass
      * False = fail
    """
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    # Find candidate dups
    idx, d2d, d3d = match_coordinates_sky(c_main, c_main, nthneighbor=2)
    cand_dups = d2d < 2*u.arcsec
    # Finish
    if np.sum(cand_dups) > 0:
        return False
    else:
        return True


def get_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table
    newdb : Table
    chk : bool, optional
      Perform some checks
    idkey : str, optional
      Key for ID

    Returns
    -------
    ids : ndarray (int)
      Old IDs are filled with negative their value
      New IDs are generated as needed

    """
    cdict = defs.get_cat_dict()
    IDs = np.zeros(len(newdb), dtype=int)
    # Setup
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    c_new = SkyCoord(ra=newdb['RA'], dec=newdb['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_new, c_main, nthneighbor=1)
    new = d2d > cdict['match_toler']
    # Old IDs
    IDs[~new] = -1 * maindb[idkey][idx[~new]]
    nnew = np.sum(new)
    # New IDs
    if nnew > 0:
        # Generate
        newID = np.max(maindb[idkey]) + 1
        newIDs = newID + np.arange(nnew, dtype=int)
        # Insert
        IDs[new] = newIDs
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    # Return
    return IDs


def set_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
    """ Set the new IDs
    Parameters
    ----------
    maindb
    newdb
    toler
    chk

    Returns
    -------
    cut_db : Table
      Cut to the new QSOs
    new : bool array
    ids : ID values

    """
    # IDs
    ids = get_new_ids(maindb, newdb, idkey=idkey)  # Includes new and old
    # Truncate
    new = ids > 0
    cut_db = newdb[new]
    cut_db.add_column(Column(ids[new], name=idkey))
    # Reset IDs to all positive
    ids = np.abs(ids)
    #
    return cut_db, new, ids


def start_maindb(private=False):
    """ Start the main DB catalog

    Returns
    -------
    maindb : Table
    tkeys : list
      List of columns in the table
    private : bool, optional
      Private DB?

    """
    idict = defs.get_db_table_format()
    #if private:
    #    idict['PRIV_ID'] = 0
        #idict.pop('IGM_ID')
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)
    # Return
    return maindb, tkeys
