#!/usr/bin/env python

""" Loads and plots a requested spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_plot script v0.3')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1")
    parser.add_argument("dbase", type=str, help="Database [igmspec,uvqs,all,priv]")
    parser.add_argument("--tol", default=5., type=float, help="Maximum offset in arcsec [default=5.]")
    #parser.add_argument("--meta", default=True, help="Show meta data? [default: True]", action="store_true")
    parser.add_argument("-g", "--group", help="Restrict on Group (e.g. BOSS_DR12)")
    parser.add_argument("-s", "--select", default=0, type=int, help="Index of spectrum to plot (when multiple exist)")
    parser.add_argument("--mplot", default=False, help="Use simple matplotlib plot [default: False]")
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
    from linetools.scripts.utils import coord_arg_to_coord

    # init
    Specdb = load_db(args.dbase, db_file=args.db_file, **kwargs)

    # Grab
    icoord = coord_arg_to_coord(args.coord)
    if args.group is not None:
        groups=[args.group]
    else:
        groups = None
    spec, meta = Specdb.spectra_from_coord(icoord, tol=args.tol*u.arcsec, groups=groups)

    # Outcome
    if meta is None:
        print("No source found, try another location or a larger tolerance.")
        return
    elif len(meta) == 1:  # One group hit
        print("Source located in group: {:s}".format(meta['GROUP'][0]))
    else:  # More than 1 spectrum
        print("More than 1 spectrum found for input source. Here is a summary:")
        indices = np.arange(len(meta)).astype(int)
        meta['INDEX'] = indices
        print(meta[['INDEX','GROUP','RA_GROUP','DEC_GROUP',Specdb.idkey,'INSTR','DISPERSER','GROUP_ID']])
        idx = args.select
        print("Plotting index={:d} which you can specify with --select".format(idx))
        #meta = all_meta[idx]
        #groups = [meta.meta['group'] for meta in all_meta]
        #print("Using group {:s}.  You can choose from this list {}".format(groups[idx], groups))

    #if args.meta:
    #    group_utils.show_group_meta(meta)

    # Load spectra
    spec.select = args.select
    if unit_test:
        return
    # Show
    if args.mplot:
        spec.plot()
    else:
        spec.plot(xspec=True)
