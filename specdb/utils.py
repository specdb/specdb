""" Module for DB utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings
import pdb

try:
    basestring
except NameError:  # For Python 3
    basestring = str

def load_db(db_type, **kwargs):
    """
    Parameters
    ----------
    db_type : str
      Type of database
      Current allowed entries are [igmspec]


    Returns
    -------
    dbobj : SpecDB object
      Generally a child

    """
    if db_type == 'igmspec':
        from specdb.specdb import IgmSpec
        Specdb = IgmSpec(**kwargs)
    elif db_type == 'priv':  # Private
        from specdb.specdb import SpecDB
        Specdb = SpecDB(**kwargs)
    else:
        raise IOError("Not ready for this dbase value: {:s}".format(db_type))

    # Return
    return Specdb


def query_table(tbl, qdict, ignore_missing_keys=True, verbose=True):
    """ Find all rows in the input table satisfying
    the query given by qdict
    Parameters
    ----------
    tbl : Table
    qdict : dict
      See query_dict documentation for rules

    Returns
    -------
    match : bool ndarray
        True = Row satisfies the query
    """
    # Init
    match = np.array([True]*len(tbl))
    tkeys = tbl.keys()

    # Perform
    for key,value in qdict.items():
        # Deal with BITWISE
        if '-BITWISE' in key:
            if '-BITWISE-OR' in key:
                flg_bitwise = 1
            elif '-BITWISE-AND' in key:
                flg_bitwise = 2
            else:
                raise IOError("Undefined BITWISE approach")
            # Truncate key name
            key = key[:key.rfind('-BITWISE')]
        else:
            flg_bitwise = 0
        # Check
        if key not in tkeys:
            msg = "Key {:s} in query_dict is not present in Table".format(key)
            if ignore_missing_keys:
                if verbose:
                    print(msg)
                continue
            else:
                raise IOError(msg)
        # Proceed
        if isinstance(value,tuple):
            # Check
            if len(value) != 2:
                raise IOError("Tuple for key={:s} in query_dict must have length 2 for min/max".format(key))
            # Min/Max
            match = match & (tbl[key].data >= value[0]) & (tbl[key].data <= value[1])
        elif isinstance(value,(list,float,basestring,int)):
            if flg_bitwise > 0: # BITWISE
                if isinstance(value,(int)):
                    match &= (tbl[key].data & 2**value).astype(bool)
                elif isinstance(value,(list)):
                    if flg_bitwise == 1:
                        bit_match = np.array([False]*len(tbl))
                    else:
                        bit_match = np.array([True]*len(tbl))
                    for item in value:
                        if flg_bitwise == 1: # OR
                            bit_match += (tbl[key].data & item).astype(bool)
                        else: # AND
                            bit_match &= (tbl[key].data & item).astype(bool)
                    match &= bit_match
            else:
                # Recast
                mlist = np.array(value).astype(tbl[key].dtype)
                # Match
                match &= np.in1d(tbl[key].data, mlist)
        else:
            raise IOError("Bad data type for query_dict value: {}".format(type(value)))
    # Return
    return match
