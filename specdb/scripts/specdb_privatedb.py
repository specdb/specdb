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
    parser.add_argument("tree_path", type=str, help="Path to the directory tree of spectral files")
    parser.add_argument("db_name", type=str, help="Name of your private DB")
    parser.add_argument("outfile", type=str, help="Filename for the private DB HDF5")
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
    from specdb.build import privatedb as pbuild

    pbuild.mk_db(pargs.tree_path, pargs.db_name, 'tmp.hdf5', ztbl, fname=True)

##
if __name__ == '__main__':
    # Giddy up
    main()
