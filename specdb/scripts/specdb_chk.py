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
    parser = argparse.ArgumentParser(description='Check a specdb DB file [v1.1]')
    parser.add_argument("db_file", type=str, help="specdb Database file (expecting an HDF5 file)")
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
        warnings.warn('DB file has no name.  Must be really old.')
        dbname = None
    else:
        print("specdb DB file is from the {:s} database".format(str(dbname)))

    # Check Creation Date
    try:
        cdate = hdf['catalog'].attrs['CREATION_DATE']
    except:
        warnings.warn('DB file has no Creation Date.  Must be really old..')
        return
    else:
        version = hdf['catalog'].attrs['VERSION']
        print("specdb DB file version={} was created on {}".format(str(version),str(cdate)))
        if dbname is not None:
            try:
                print("Latest version for specdb DB type={} is version={}".format(
                        dbname, str(dbinfo[dbname]['latest_version'])))
            except KeyError:
                print("No version version history for {:s}".format(str(dbname)))
            else:
                # Check
                print("Latest creation date for this DB version was {}".format(
                    dbinfo[dbname][version]['newest_date']))
                print("Oldest valid DB file for this DB version was {}".format(
                        dbinfo[dbname][version]['oldest_ok_date']))
                # Compare?

    # Sources and spectra
    print("-------------------------------------------------")
    print("There are {:d} unique sources in the source catalog".format(len(hdf['catalog'].value)))

    # List datasets
    nspec = 0
    for key in hdf.keys():
        if key == 'catalog':
            continue
        print("Dataset: {:s}".format(key))
        # Spectra
        try:
            nspec += hdf[key]['spec'].size
        except:
            pass

    print("There are a total of {:d} spectra in the database".format(nspec))

##
if __name__ == '__main__':
    # Giddy up
    main()
