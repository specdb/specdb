""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky


def zem_from_radec(ra, dec, catalog, toler=2*u.arcsec, debug=False):
    """ Parse input catalog (e.g. Myers) for zem

    Parameters
    ----------
    ra : list or array
      RA in deg
    dec : list or array
      DEC in deg
    catalog : Table
      Must contain RA,DEC,ZEM,ZEM_SOURCE
    debug : bool, optional

    Returns
    -------
    zem : array
      Redshifts
    zsource : array
      str array of sources
    """
    # Generate coordinates
    icoord = SkyCoord(ra=ra, dec=dec, unit='deg')
    # Quasar catalog
    qcoord = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], unit='deg')
    # Match
    idx, d2d, d3d = match_coordinates_sky(icoord, qcoord, nthneighbor=1)
    good = d2d < toler
    if debug:
        pdb.set_trace()
    # Finish
    zem = np.zeros(len(ra))
    try:
        zem[good] = catalog['ZEM'][idx[good]]
    except IndexError:
        pdb.set_trace()
    zsource = np.array([str('NONENONE')]*len(ra))
    zsource[good] = catalog['ZEM_SOURCE'][idx[good]]

    # Return
    return zem, zsource

