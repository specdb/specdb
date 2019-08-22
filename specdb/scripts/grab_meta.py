#!/usr/bin/env python

""" Grab and show all the meta info for a source
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_meta script v0.1')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1")
    parser.add_argument("dbase", type=str, help="Database [igmspec,uvqs,qpq,casbah,priv]")
    parser.add_argument("--tol", default=5., type=float, help="Maximum offset in arcsec [default=5.]")
    parser.add_argument("-g", "--group", help="Restrict on Group (e.g. BOSS_DR12)")
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
    if args.group is not None:
        groups=[args.group]
    else:
        groups = None
    meta = Specdb.meta_from_position(icoord, args.tol*u.arcsec, groups=groups)
    if unit_test:
        return meta

    # Outcome
    if meta is None:
        print("No source found within the data groups, try another location or a larger tolerance.")
        return
    else:
        group_utils.show_group_meta(meta, idkey=Specdb.idkey, show_all_keys=False)

