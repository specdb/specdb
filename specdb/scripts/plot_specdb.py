#!/usr/bin/env python

""" Loads and plots a requested spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_plot script v0.3')
    parser.add_argument("coord", type=str, help="Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322 or 07:45:00.47,34:17:31.1")
    parser.add_argument("dbase", type=str, help="Database [igmspec,all,priv]")
    parser.add_argument("--tol", default=5., type=float, help="Maximum offset in arcsec [default=5.]")
    parser.add_argument("--meta", default=True, help="Show meta data? [default: True]", action="store_true")
    parser.add_argument("-s", "--survey", help="Name of Survey to use")
    parser.add_argument("--select", default=0, type=int, help="Index of spectrum to plot (when multiple exist)")
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

    from astropy import units as u
    from specdb.utils import load_db
    from linetools.scripts.utils import coord_arg_to_coord

    # init
    Specdb = load_db(args.dbase, db_file=args.db_file, **kwargs)

    # Grab
    icoord = coord_arg_to_coord(args.coord)
    all_spec, all_meta = Specdb.spec_from_coord(icoord, tol=args.tol*u.arcsec, isurvey=args.survey)

    # Outcome
    if len(all_meta) == 0:
        print("No source found, try another location or a larger tolerance.")
        return
    elif len(all_meta) == 1:  # One survey hit
        spec = all_spec[0]
        meta = all_spec[0]
    else:  # More than 1 survey
        idx = 0
        spec = all_spec[idx]
        meta = all_meta[idx]
        surveys = [meta.meta['survey'] for meta in all_meta]
        print("Source located in more than one survey")
        print("Using survey {:s}.  You can choose from this list {}".format(surveys[idx], surveys))

        #print("Choose another survey from this list (as you wish): {}".format(surveys))

    if args.meta:
        Specdb.idb.show_meta(imeta=meta)

    # Load spectra
    spec.select = args.select
    if unit_test:
        return
    # Show  [may transition to xspec]
    if args.mplot:
        spec.plot()
    else:
        spec.plot(xspec=True)
