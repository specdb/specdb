#!/usr/bin/env python

""" Loads (and can plot) a requested SDSS/BOSS spectrum
Default is by PLATE/FIBER
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='specdb_sdss script v0.1')
    parser.add_argument("plate", type=int, help="Plate")
    parser.add_argument("fiberid", type=int, help="FiberID")
    parser.add_argument("dbase", type=str, help="Database [igmspec,all]")
    parser.add_argument("-s", "--survey", help="Name of Survey to use (BOSS_DR12 or SDSS_DR7)")
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
    Specdb = load_db(args.dbase, **kwargs)

    # Load meta table(s)
    if args.survey is None:
        surveys = ['BOSS_DR12', 'SDSS_DR7']
    else:
        surveys = [args.survey]
    for kk,survey in enumerate(surveys):
        meta = Specdb[survey].meta
        if 'FIBERID' not in meta.keys():
            meta.rename_column('FIBER','FIBERID')
        if kk > 0:
            mtbl = vstack([mtbl, meta], join_type='inner')
        else:
            mtbl = meta

    # Find plate/fiber
    imt = np.where((mtbl['PLATE'] == args.plate) & (mtbl['FIBERID'] == args.fiberid))[0]
    if len(imt) == 0:
        print("Plate and Fiber not found.  Try again")
        return
    else:
        mt = imt[0]
        scoord = SkyCoord(ra=mtbl['RA_GROUP'][mt], dec=mtbl['DEC_GROUP'][mt], unit='deg')

    # Grab
    print("Grabbing data for J{:s}{:s}".format(scoord.ra.to_string(unit=u.hour,sep='',pad=True),
                                              scoord.dec.to_string(sep='',pad=True,alwayssign=True)))
    spec, meta = Specdb.spectra_from_coord(scoord, groups=surveys)

    # Outcome
    if len(meta) == 0:
        print("No source found, try another location or a larger tolerance.")
        return
    elif len(meta) == 1:  # One survey hit
        pass
    else:  # More than 1 survey
        pdb.set_trace()
        idx = 0
        spec = all_spec[idx]
        meta = all_meta[idx]
        surveys = [meta.meta['group'] for meta in all_meta]
        print("Source located in more than one SDSS survey")
        print("Using survey {:s}".format(surveys[idx]))
        print("You can choose from this list {}".format(surveys))

    #group_utils.show_group_meta(meta)

    # Load spectra
    spec.select = args.select
    if unit_test:
        return
    # Show  [may transition to xspec]
    if args.plot:
        spec.plot(xspec=True)
