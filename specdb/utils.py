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


def clean_vstack(tables, labels, **kwargs):
    """ Perform an astropy.table.vstack on a list of Tables
    after first renaming any conflicting columns

    Parameters
    ----------
    tables : list
     list of astropy.table.Table
    labels : list
     list of str labels to append to unruly column names
    kwargs

    Returns
    -------
    stack : Table

    """
    from astropy.table import operations, TableMergeError, vstack

    # Loop on tables
    tkeys = []
    tbl_idx = []
    for ss, table in enumerate(tables):
        for key in table.keys():
            try:
                idx = tkeys.index(key)
            except ValueError:
                tkeys.append(key)
                tbl_idx.append(ss)
            else:
                rename = False
                # Check data type
                try:
                    operations.common_dtype([table[key],tables[tbl_idx[idx]][key]])
                except TableMergeError:
                    rename = True
                # Check shape
                shapes = [table[key].shape[1:],tables[tbl_idx[idx]][key].shape[1:]]
                uniq_shapes = set(shapes)
                if len(uniq_shapes) != 1:
                    rename = True
                # Rename?
                if rename:
                    new_name = labels[ss]+'_'+key
                    warnings.warn("Renaming column in table {:s} from {:s} to {:s}".format(labels[ss], key, new_name))
                    table.rename_column(key, new_name)
                    tkeys.append(new_name)
                    tbl_idx.append(ss)
    # Stack
    return vstack(tables)

def hdf_decode(obj, itype=None):
    """ Decode the incoming hdf5 object
    Usually byte to str
    Parameters
    ----------
    obj : object

    Returns
    -------
    dobj : object
      decoded object

    """
    if itype == 'Table':
        from astropy.table import Table, Column
        dobj = Table(obj)
        # FIX STRING COLUMNS
        for key in dobj.keys():
            if 'bytes' in dobj[key].dtype.name:
                ss = [hdf_decode(ii) for ii in dobj[key]]
                # Remove
                dobj.remove_column(key)
                # Add back
                dobj[key] = Column(ss)
    else:  # Auto
        if isinstance(obj,bytes):
            dobj = obj.decode('utf-8')
        else:
            dobj = obj
        #
    return dobj

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
    elif db_type == 'uvqs':
        from specdb.specdb import UVQS
        Specdb = UVQS(**kwargs)
    elif db_type == 'priv':  # Private
        from specdb.specdb import SpecDB
        Specdb = SpecDB(**kwargs)
    else:
        raise IOError("Not ready for this dbase value: {:s}".format(db_type))

    # Return
    return Specdb


def query_table(tbl, qdict, ignore_missing_keys=True, verbose=True,
                tbl_name=''):
    """ Find all rows in the input table satisfying
    the query given by qdict
    Parameters
    ----------
    tbl : Table
    qdict : dict
      See query_dict documentation for rules
    ignore_missing_keys : bool, optional
      Ignores any keys in the query_dict not found in Table
      Otherwise, throw an IOError
    tbl_name : str, optional
      Name of table.  Mainly for error message

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
            msg = "Key {:s} in query_dict is not present in Table {:s}".format(key, tbl_name)
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
                    try:
                        match &= bit_match
                    except TypeError:
                        pdb.set_trace()
            else:
                # Recast
                mlist = np.array(value).astype(tbl[key].dtype)
                # Match
                match &= np.in1d(tbl[key].data, mlist)
        else:
            raise IOError("Bad data type for query_dict value: {}".format(type(value)))
    # Return
    return match
