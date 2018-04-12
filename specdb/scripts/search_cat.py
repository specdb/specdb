#!/usr/bin/env python

""" Grab and show all the catalog info for a source
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_cat script v0.1')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1")
    parser.add_argument("dbase", type=str, help="Database [igmspec,uvqs,qpq,priv]")
    parser.add_argument("--tol", default=5., type=float, help="Maximum offset in arcsec [default=5.]")
    parser.add_argument("--db_file", help="Full path of db_file")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np
    from astropy import units as u
    from specdb.utils import load_db
    from specdb import group_utils
    from linetools.scripts.utils import coord_arg_to_coord

    # init
    Specdb = load_db(args.dbase, db_file=args.db_file, **kwargs)


    # Grab
    icoord = coord_arg_to_coord(args.coord)

    # Cut down using source catalog
    matches, sub_cat, IDs = Specdb.qcat.query_position(icoord, args.tol*u.arcsec)

    print(sub_cat['RA','DEC','zem','flag_group','STYPE'])

