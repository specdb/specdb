#!/usr/bin/env python
"""
Check a DB file from specdb
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Check a specdb DB file')
    parser.add_argument("db_file", type=str, help="Database file")
    #parser.add_argument("-v", "--version", default='v01', help="DB version to generate")
    #parser.add_argument("-v", "--version", default='v01', help="DB version to generate")
    #parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    Parameters
    ----------
    pargs

    Returns
    -------

    """
    import h5py
    import warnings
    from specdb import defs
    dbinfo = defs.dbase_info()

    # Open file
    hdf = h5py.File(pargs.db_file, 'r')

    # Name
    try:
        dbname = hdf['catalog'].attrs['NAME']
    except:
        warnings.warn('DB file has no name.  Must be really old..')
        dbname = None
    else:
        print("specdb DB file is from the {:s} database".format(dbname))

    # Check Creation Date
    try:
        cdate = hdf['catalog'].attrs['CREATION_DATE']
    except:
        warnings.warn('DB file has no Creation Date.  Must be really old..')
        return
    else:
        version = hdf['catalog'].attrs['VERSION']
        print("specdb DB file version={:s} was created on {:s}".format(version,cdate))
        if dbname is not None:
            try:
                print("Latest version for specdb DB type={:s} is version={:s}".format(
                        dbname, dbinfo[dbname]['latest_version']))
            except KeyError:
                print("No version version history for {:s}".format(dbname))
            else:
                # Check
                print("Latest creation date for this DB version was {:s}".format(
                    dbinfo[dbname][version]['newest_date']))
                print("Oldest valid DB file for this DB version was {:s}".format(
                        dbinfo[dbname][version]['oldest_ok_date']))
                # Compare?

    # List datasets
    for key in hdf.keys():
        print("Dataset: {:s}".format(key))

##
if __name__ == '__main__':
    # Giddy up
    main()
