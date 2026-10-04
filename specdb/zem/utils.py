""" Module for catalog utilities
"""

import numpy as np

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
    catalog : Table or None
      Must contain RA,DEC,ZEM,ZEM_SOURCE
    debug : bool, optional

    Returns
    -------
    zem : array
      Redshifts
    zsource : array
      str array of sources
    ZQ : array, optional
    """
    # For initializing without redshifts
    if catalog is None:
        nobj = len(ra)
        return -1*np.ones_like(ra), np.array(['UNKN']*nobj), np.zeros(nobj)

    # Generate coordinates
    icoord = SkyCoord(ra=ra, dec=dec, unit='deg')
    # Quasar catalog
    qcoord = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], unit='deg')
    # Match
    idx, d2d, d3d = match_coordinates_sky(icoord, qcoord, nthneighbor=1)
    good = d2d < toler
    # Finish
    zem = np.zeros(len(ra))
    try:
        zem[good] = catalog['ZEM'][idx[good]]
    except IndexError:
        raise IndexError(
            "Redshift catalog is too short for the matched indices; "
            "check that 'catalog' contains RA, DEC, ZEM and ZEM_SOURCE.")
    # An explicit width: np.array(['NONENONE']*n) is '<U8', which silently
    # truncated any longer source name -- 'BOSS_PIPE' (a real z_priority entry,
    # 9 characters) became 'BOSS_PIP'.  dtype=object would fix that too but
    # h5py cannot write object arrays.
    zsource = np.array(['NONENONE']*len(ra), dtype='<U32')
    zsource[good] = catalog['ZEM_SOURCE'][idx[good]]

    # ZQ?
    if 'ZQ' in catalog.keys():
        ZQ = np.zeros_like(zem, dtype=int)
        ZQ[good] = catalog['ZQ'][idx[good]]
        return zem, zsource, ZQ
    else:
        return zem, zsource, None

