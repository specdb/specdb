""" Module of utilities for building DB files
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings

import numbers
import pdb

from igmspec import defs

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

#from linetools import utils as ltu



def add_to_flag(cur_flag, add_flag):
    """ Add a bitwise flag to an existing flat

    Parameters
    ----------
    cur_flag : int or ndarray
    add_flag : int

    Returns
    -------
    new_flag : int or ndarray

    """
    if isinstance(cur_flag, numbers.Integral):
        if (cur_flag % add_flag) < (add_flag//2):
            cur_flag += add_flag
        return cur_flag
    else:  # Array
        mods = cur_flag % add_flag
        upd = mods < (add_flag//2)
        cur_flag[upd] += add_flag
        return cur_flag


def chk_maindb_join(maindb, newdb):
    """Check that new data is consistent with existing table

    Parameters
    ----------
    maindb : Table
    newdb : Table

    Returns
    -------
    answer : bool

    """
    # One way
    for key in newdb.keys():
        try:
            assert key in maindb.keys()
        except AssertionError:
            pdb.set_trace()
            return False
    return True


def chk_for_duplicates(maindb):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table

    Return
    ------
    result : bool
      * True = pass
      * False = fail
    """
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    # Find candidate dups
    idx, d2d, d3d = match_coordinates_sky(c_main, c_main, nthneighbor=2)
    cand_dups = d2d < 2*u.arcsec
    # Finish
    if np.sum(cand_dups) > 0:
        return False
    else:
        return True


def get_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table
    newdb : Table
    chk : bool, optional
      Perform some checks
    idkey : str, optional
      Key for ID

    Returns
    -------
    ids : ndarray (int)
      Old IDs are filled with negative their value
      New IDs are generated as needed

    """
    cdict = defs.get_cat_dict()
    IDs = np.zeros(len(newdb), dtype=int)
    # Setup
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    c_new = SkyCoord(ra=newdb['RA'], dec=newdb['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_new, c_main, nthneighbor=1)
    new = d2d > cdict['match_toler']
    # Old IDs
    IDs[~new] = -1 * maindb[idkey][idx[~new]]
    nnew = np.sum(new)
    # New IDs
    if nnew > 0:
        # Generate
        newID = np.max(maindb[idkey]) + 1
        newIDs = newID + np.arange(nnew, dtype=int)
        # Insert
        IDs[new] = newIDs
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    # Return
    return IDs


def set_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
    """ Set the new IDs
    Parameters
    ----------
    maindb
    newdb
    toler
    chk

    Returns
    -------
    cut_db : Table
      Cut to the new QSOs
    new : bool array
    ids : ID values

    """
    # IDs
    ids = get_new_ids(maindb, newdb, idkey=idkey)  # Includes new and old
    # Truncate
    new = ids > 0
    cut_db = newdb[new]
    cut_db.add_column(Column(ids[new], name=idkey))
    # Reset IDs to all positive
    ids = np.abs(ids)
    #
    return cut_db, new, ids


def start_maindb(private=False):
    """ Start the main DB catalog

    Returns
    -------
    maindb : Table
    tkeys : list
      List of columns in the table
    private : bool, optional
      Private DB?

    """
    idict = defs.get_db_table_format()
    #if private:
    #    idict['PRIV_ID'] = 0
        #idict.pop('IGM_ID')
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)
    # Return
    return maindb, tkeys


def chk_meta(meta):#, skip_igmid=False):
    """ Vettes a meta Table prior to its being ingested into the hdf

    Parameters
    ----------
    meta

    Returns
    -------
    chk : bool

    """
    from specdb.defs import instruments, get_req_clms
    from astropy.time import Time
    from astropy.table import Column
    # Init
    inst_dict = instruments()

    chk = True
    # Required columns
    req_clms = get_req_clms()
    meta_keys = meta.keys()
    for clm in req_clms:
        if clm not in meta_keys:
            chk = False
            print("Missing column {:s} in meta".format(clm))
    # Check date formatting
    try:
        tval = Time(list(meta['DATE-OBS'].data), format='iso')
    except:
        print("Bad DATE-OBS formatting")
        chk = False
    # Check for unicode
    for key in meta_keys:
        if 'unicode' in meta[key].dtype.name:
            warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))
            tmp = Column(meta[key].data.astype(str), name=key)
            meta.remove_column(key)
            meta[key] = tmp
    # Check instrument
    meta_instr = meta['INSTR'].data
    db_instr = np.array(inst_dict.keys()).astype(str)
    if not np.all(np.in1d(meta_instr, db_instr)):
        print("Bad instrument in meta data")
        chk = False
    # Return
    return chk


def set_resolution(head, instr=None):
    """ Sets resolution based on the instrument and header

    Parameters
    ----------
    head : FITS header
    instr : str, optional
      If not provided, attempt to grab from header

    Returns
    -------

    """
    from igmspec import defs
    # Dicts
    Rdicts = defs.get_res_dicts()
    # Grab instrument
    if instr is None:
        if 'CURRINST' in head.keys():  # ESI, NIRSPEC
            instr = head['CURRINST'].strip()
        elif 'INSTRUME' in head.keys():
            if 'HIRES' in head['INSTRUME']:
                instr = 'HIRES'
            elif 'MOSFIRE' in head['INSTRUME']:
                instr = 'MOSFIRE'
            elif 'MIKE' in head['INSTRUME']:
                instr = 'MIKE'
            elif 'MagE' in head['INSTRUME']:
                instr = 'MagE'
            elif 'GMOS' in head['INSTRUME']:
                instr = 'GMOS'
            elif 'GNIRS' in head['INSTRUME']:
                instr = 'GNIRS'
            elif 'NIRI' in head['INSTRUME']:
                instr = 'NIRI'
            elif 'mmt' in head['INSTRUME']:
                instr = 'mmt'
            elif 'MODS1B' in head['INSTRUME']:
                instr = 'MODS1B'
            elif 'MODS1R' in head['INSTRUME']:
                instr = 'MODS1R'
        else:
            pass
        if instr is None:
            raise ValueError("NEED MORE INFO FOR INSTR")

    # Grab resolution
    if instr == 'ESI':
        try:
            return Rdicts[instr][head['SLMSKNAM']]
        except KeyError:
            pdb.set_trace()
    elif instr == 'HIRES':
        try:
            return Rdicts[instr][head['DECKNAME'].strip()]
        except KeyError:
            print("Need to add {:s}".format(head['DECKNAME']))
            pdb.set_trace()
    elif instr == 'GMOS':
        try:
            return Rdicts[instr][head['GRATING']]
        except KeyError:
            print("Need to add {:s}".format(head['GRATING']))
            pdb.set_trace()
    elif instr == 'MOSFIRE':
        try:
            res = Rdicts[instr][head['FILTER']]*0.7
        except KeyError:
            print("Need to add {:s}".format(head['FILTER']))
            pdb.set_trace()
        else:
            swidth = defs.slit_width(head['MASKNAME'])
            return res/swidth
    elif instr == 'GNIRS':
        try:
            res = Rdicts[instr][head['GRATING']]*0.3
        except KeyError:
            print("Need to add {:s}".format(head['GRATING']))
            pdb.set_trace()
        else:
            swidth = defs.slit_width(head['SLIT'])
            return res/swidth
    elif instr == 'NIRI':
        try:
            res = Rdicts[instr][head['FILTER3']]/4.
        except KeyError:
            print("Need to add {:s}".format(head['FILTER3']))
            pdb.set_trace()
        else:
            swidth = defs.slit_width(head['FPMASK'])  #PIXELS
            return res*swidth
    elif instr == 'NIRSPEC':  # LOW DISPERSION
        try:
            return 2000.*0.38/defs.slit_width(head['SLITNAME'])
        except KeyError:
            print("Need to add {:s}".format(head['SLITNAME']))
            pdb.set_trace()
    elif instr == 'MagE':
        try:
            return 4100./defs.slit_width(head['SLITNAME'])
        except KeyError:
            print("Need to add {:s}".format(head['SLITNAME']))
            pdb.set_trace()
    elif 'mmt' in instr:
        try:
            res = Rdicts[instr][head['DISPERSE']]*0.6
        except KeyError:
            print("Need to add {:s}".format(head['DISPERSE']))
            pdb.set_trace()
        else:
            swidth = defs.slit_width(head['APERTURE'])
            return res/swidth
    elif instr in ['MODS1B','MODS1R']:
        try:
            res = Rdicts[instr][head['GRATNAME']]*0.6
        except KeyError:
            print("Need to add {:s}".format(head['GRATNAME']))
            pdb.set_trace()
        else:
            swidth = defs.slit_width(head['MASKNAME'])
            return res/swidth
    elif instr == 'MIKE':
        try:
            res = Rdicts[head['INSTRUME']]
        except KeyError:
            print("Need to add {:s}".format(instr))
            pdb.set_trace()
        else:
            try:
                swidth = defs.slit_width(head['SLITSIZE'])
            except TypeError:
                warnings.warn("MIKE slit not given in header. Assuming 1.0")
                swidth = 1.
            return res/swidth
    else:
        raise IOError("Not read for this instrument")
