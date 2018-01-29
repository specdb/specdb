#!/usr/bin/env python

""" Loads (and can plot) a requested SDSS/BOSS spectrum
Default is by PLATE/FIBER
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_sdss script v0.2')
    parser.add_argument("plate", type=int, help="Plate")
    parser.add_argument("fiberid", type=int, help="FiberID")
    parser.add_argument("--dbase", default='igmspec', type=str, help="Database [igmspec,all]")
    parser.add_argument("--survey", help="Name of Survey to use (BOSS_DR12 or SDSS_DR7)")
    parser.add_argument("--select", default=0, type=int, help="Index of spectrum to plot (when multiple exist)")
    parser.add_argument("-p", "--plot", default=False, action="store_true", help="Plot with lt_xspec")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, **kwargs):
    """ Run
    """
    import numpy as np
    from astropy.table import vstack
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    from specdb.utils import load_db
    from specdb import group_utils

    # init
    if args.dbase not in ['igmspec']:
        raise IOError("Not prepared for this database: {:s}".format(args.dbase))
    Specdb = load_db(args.dbase, **kwargs)

    if args.survey is None:
        surveys = ['BOSS_DR12', 'SDSS_DR7']
    else:
        surveys = [args.survey]

    # Load
    spec, meta = Specdb.get_sdss(args.plate, args.fiberid, groups=surveys)
    if spec is None:
        if meta == -1:
            return

    # Outcome
    if len(meta) == 0:
        print("No source found, try another location or a larger tolerance.")
        return
    elif len(meta) == 1:  # One survey hit
        pass
    else:  # More than 1 survey
        surveys = [meta['GROUP'] for row in meta]
        print("Source located in more than one SDSS survey")
    indices = np.arange(len(meta)).astype(int)
    meta['INDEX'] = indices
    print(meta[['INDEX','GROUP','RA_GROUP','DEC_GROUP',Specdb.idkey,'INSTR','DISPERSER','GROUP_ID']])

    # Add labels
    lbls = []
    for imeta in meta:
        lbls.append('{:d}_{:s}'.format(imeta['INDEX'], imeta['GROUP']))
    spec.labels = lbls
    # Load spectra
    spec.select = args.select
    if unit_test:
        return
    # Show  [may transition to xspec]
    if args.plot:
        spec.plot(xspec=True)
