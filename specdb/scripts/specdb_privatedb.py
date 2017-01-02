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
    parser.add_argument("--version", type=str, help="Version of the DB; default is `v00`")
    parser.add_argument("--publisher", type=str, help="Publisher of the DB; default is `Unknown`")
    parser.add_argument("--fname", default=False, action="store_true", help="Parse RA/DEC from filename?")

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
        raise NotImplementedError("Not ready for this yet.  You should add your table to the database tree with extension _ztbl.fits")
    if pargs.zspecdb not in [None, 'igmspec']:
        raise NotImplementedError("Not ready for this yet")

    # define main tree
    tree = pargs.tree_path

    # iztbl
    if pargs.zspecdb == 'igmspec':
        iztbl = 'igmspec'
    else:
        # Search for a z table
        ztbl_files = glob.glob(tree+'/*_ztbl*')
        if len(ztbl_files) == 1:
            print("Reading redshift table {:s}".format(ztbl_files[0]))
            iztbl = Table.read(ztbl_files[0])
        elif len(ztbl_files) == 0:
            raise IOError("No redshift table provided")
        else:
            raise IOError("Multiple redshift tables found in your tree")

    # version
    if pargs.version is None:
        version = 'v00'
    else:
        version = pargs.version
    # version
    if pargs.publisher is None:
        publisher = 'Unknown'
    else:
        publisher = pargs.publisher

    # Run
    pbuild.mk_db(pargs.db_name, tree, pargs.outfile, iztbl,
                 fname=pargs.fname, version=version, publisher=publisher)

##
if __name__ == '__main__':
    # Giddy up
    main()
