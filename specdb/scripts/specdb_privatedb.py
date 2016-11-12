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
    parser = argparse.ArgumentParser(description='Generate a private specdb DB')
    parser.add_argument("db_name", type=str, help="Name of your private DB")
    parser.add_argument("tree_path", type=str, help="Path to the directory tree of spectral files")
    parser.add_argument("outfile", type=str, help="Filename for the private DB HDF5")
    parser.add_argument("--ztbl", help="Name of data file containing redshift info")
    parser.add_argument("--zspecdb", help="Name of specdb DB to use for redshifts")

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
    import glob
    from astropy.table import Table
    from specdb.build import privatedb as pbuild

    if pargs.ztbl is not None:
        raise NotImplementedError("Not ready for this yet")
    if pargs.zspecdb is not None:
        raise NotImplementedError("Not ready for this yet")

    # Search for a z table
    tree = pargs.tree_path
    ztbl_files = glob.glob(tree+'/*_ztbl*')
    if len(ztbl_files) == 1:
        print("Reading redshift table {:s}".format(ztbl_files[0]))
        ztbl = Table.read(ztbl_files[0])
    elif len(ztbl_files) == 0:
        raise IOError("No redshift table provided")
    else:
        raise IOError("Multiple redshift tables found in your tree")
    # Run
    pbuild.mk_db(pargs.db_name, tree, pargs.outfile, ztbl, fname=True)

##
if __name__ == '__main__':
    # Giddy up
    main()
