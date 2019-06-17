""" Module of utilities for building DB files
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import warnings

import numbers
import pdb

from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.time import Time

from linetools import utils as ltu

from specdb import defs
from specdb.cat_utils import match_ids
from specdb.utils import clean_vstack

try:
    bstr = bytes
except NameError:  # For Python 2
    bstr = str

def add_ids(maindb, meta, flag_g, tkeys, idkey, first=False, debug=False, **kwargs):
    """ Add IDs to
    Input meta table has its CAT_ID values set in place

    Parameters
    ----------
    maindb : Table
      Main catalog
    meta : Table
      Meta table being added
    flag_g : int
      Flag for the new group
    tkeys : list
      List of main keys for the catalog
    ikdey : str
      ID key
    first : bool, optional
      First call to the routine?

    Returns
    -------
    maindb : Table
      Updated catalog table

    """
    newcut, new, ids = set_new_ids(maindb, meta, idkey, first=first, debug=debug, **kwargs)
    # If new sources
    if np.sum(new) > 0:
        newcut['flag_group'] = flag_g
        newcut.rename_column('RA_GROUP', 'RA')
        newcut.rename_column('DEC_GROUP', 'DEC')
        newcut.rename_column('zem_GROUP', 'zem')
        cat_meta = newcut[tkeys]
    else:
        print("No new sources in this group")
    # Set or append
    if first:
        maindb = cat_meta
    else:
        # Update group flags
        old_ids = ids[~new]
        midx = match_ids(old_ids, maindb[idkey].data) # np.array(maindb[idkey][ids[~new]])
        maindb['flag_group'][midx] += flag_g   # ASSUMES NOT SET ALREADY
        if np.sum(new) > 0:
            # Catalog
            assert chk_maindb_join(maindb, cat_meta)
            # Append
            maindb = vstack([maindb,cat_meta], join_type='exact')
    # Return
    return maindb


def add_to_flag(cur_flag, add_flag):
    """ Add a bitwise flag to an existing flag

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


def add_to_group_dict(group_name, gdict, skip_for_debug=False):
    """ Add input group_name to the group dict
    Done in place

    Parameters
    ----------
    group_name : str
    gdict : dict
      Dict of data groups

    Returns
    -------
    new_flag : int
      New bitwise flag
    """
    if group_name in gdict.keys():
        if not skip_for_debug:
            raise IOError("Group already exists in group dict.  Should not be here..")
    # Find new flag
    if len(list(gdict.keys())) == 0: # First one
        new_flag = 1
    else:
        max_flag = max(gdict.values())
        new_flag = 2*max_flag
    # Insert
    gdict[group_name] = new_flag
    return new_flag


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


def chk_for_duplicates(maindb, tol=2*u.arcsec, dup_lim=0):
    """ Check for duplicates in the catalog to within tol

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
    cand_dups = d2d < tol
    # Finish
    if np.sum(cand_dups) > dup_lim:
        return False
    else:
        return True


def chk_meta(meta, chk_cat_only=False):
    """ Vettes a meta Table prior to its being ingested into the hdf

    Parameters
    ----------
    meta
    chk_cat_only : bool, optional
      Check that the primary catalog keys are present
      These should match well with the ones in sdbbu.start_maindb()

    Returns
    -------
    chk : bool

    """
    # Init
    inst_dict = defs.instruments()

    chk = True
    if chk_cat_only: # Only check for catalog generation
        cat_dict = defs.get_db_table_format()
        for key in cat_dict.keys():
            if key in ['flag_group']:
                continue
            # RA/DEC are special
            if key in ['RA', 'DEC', 'zem']:
                akey = key+'_GROUP'
            else:
                akey = key
            if akey not in meta.keys():
                print("key={:s} not in meta table!".format(akey))
                chk = False
    else:
        # Required columns for main meta
        req_clms = defs.get_req_clms()
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
            pdb.set_trace()
            chk = False
        # Clean
        clean_table_for_hdf(meta)
        # Check instrument
        meta_instr = meta['INSTR'].data
        db_instr = np.array(list(inst_dict.keys())).astype(bstr)  # bytes
        if not np.all(np.in1d(meta_instr, db_instr)):
            print("Bad instrument in meta data")
            chk = False
        # ID key
        if sv_idkey not in meta_keys:
            print("Missing ID key: {:s}".format(sv_idkey))
            chk = False
    # Return
    return chk


def chk_vstack(hdf):
    """ Check whether the meta data in a specdb database can
    be stacked with specdb.utils.clean_vstack

    Parameters
    ----------
    hdf : HDF5 pointer

    Returns
    -------
    chk : bool

    """
    meta_tables = []
    labels = []
    for key in hdf.keys():
        try:
            meta = Table(hdf[key]['meta'].value)
        except (KeyError,ValueError):
            print("Skipping data group {:s}".format(key))
        else:
            # Save a snippet
            meta_tables.append(meta[0:1])
            labels.append(key)
    # Try to stack
    try:
        stack = clean_vstack(meta_tables, labels)
    except:
        chk = False
    else:
        print("Passing chk_vstack...")
        chk = True
    # Return
    return chk


def clean_table_for_hdf(tbl):
    """ Prepare an input table for writing in an HDF5 object
    Cleans unicode

    Parameters
    ----------
    tbl : Table

    Returns
    -------
    tbl is modified in place

    """
    #
    tbl_keys = tbl.keys()
    for key in tbl_keys:
        # Check for unicode
        #if 'unicode' in tbl[key].dtype.name:
        if 'U' in tbl[key].dtype.__repr__():
            warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))
            tmp = Column(tbl[key].data.astype(bstr), name=key)
            tbl.remove_column(key)
            tbl[key] = tmp

def get_new_ids(maindb, newdb, idkey, chk=True, mtch_toler=None, pair_sep=0.5*u.arcsec,
                close_pairs=False, debug=False):
    """ Generate new CAT_IDs for an input DB

    Parameters
    ----------
    maindb : Table
    newdb : Table
      RA, DEC assumed to be given by RA_GROUP and DEC_GROUP
    chk : bool, optional
      Perform some checks
    idkey : str
      Key for ID
    mtch_toler : Quantity, optional
      Matching tolerance;  typically taken from the default
    pair_sep : Angle, optional
      Sepration at which a pair is considered 'real'
    close_pairs : bool, optional
      Input list includes close pairs (i.e. within mtch_toler)

    Returns
    -------
    ids : ndarray (int)
      Old IDs are filled with negative their value
      New IDs are generated as needed

    """
    if mtch_toler is None:
        cdict = defs.get_cat_dict()
        mtch_toler = cdict['match_toler']
    IDs = np.zeros(len(newdb), dtype=int)
    # Setup
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    c_new = SkyCoord(ra=newdb['RA_GROUP'], dec=newdb['DEC_GROUP'], unit='deg')
    # Check for pairs in the new list
    pidx1, pidx2, pd2d, _ = c_new.search_around_sky(c_new, mtch_toler)
    pairs = pd2d > pair_sep
    if np.sum(pairs) and (not close_pairs):
        print ("Input catalog includes pairs closer than {:g} and wider than {:g}".format(mtch_toler, pair_sep))
        raise IOError("Use close_pairs=True if appropriate")
    # Find new sources (ignoring pairs at first)
    idx, d2d, d3d = match_coordinates_sky(c_new, c_main, nthneighbor=1)
    new = d2d > mtch_toler
    # Old IDs
    IDs[~new] = -1 * maindb[idkey][idx[~new]]
    # Now deal with pairs
    if np.sum(pairs) > 0:
        # Check against catalog
        pidx, pd2d, _ = match_coordinates_sky(c_new[pidx1][pairs], c_main, nthneighbor=1)
        not_pair_match = pd2d > pair_sep
        # Reset new -- It will get a new ID below -- np.where is needed to actually set new
        new[pidx1[pairs][np.where(not_pair_match)[0]]] = True
    # New IDs
    nnew = np.sum(new)
    new_idx = np.where(new)[0]
    newID = np.max(maindb[idkey])
    # Ingest
    if nnew == 1:
        IDs[new_idx] = newID + 1
    elif nnew > 1: # Deal with duplicates
        sub_c_new = c_new[new]
        dup_idx, dup_d2d, _ = match_coordinates_sky(sub_c_new, sub_c_new, nthneighbor=2)
        if close_pairs:
            dups = dup_d2d < pair_sep
        else:
            dups = dup_d2d < mtch_toler
        ndups = np.sum(dups)
        # Not duplicates
        IDs[new_idx[~dups]] = newID + 1 + np.arange(np.sum(~dups))
        newID = max(np.max(IDs), newID)

        # Duplicates
        if ndups > 0:
            warnings.warn("We found {:d} duplicates (e.g. multiple spectra). Hope this was expected".format(ndups//2))
            # Cut down to unique and restrict to new ones (there are at least 2 duplicates per match)
            dup_idx = np.where(dups)[0]
            dup_filled = np.array([False]*len(sub_c_new))
            if debug:
                pdb.set_trace()
            for idup in dup_idx: # Ugly loop..
                if dup_filled[idup]:  # Already filled as a duplicate
                    continue
                dcoord = sub_c_new[idup]
                sep = dcoord.separation(sub_c_new)
                if close_pairs:
                    isep = np.where(sep < pair_sep)[0]
                else:
                    isep = np.where(sep < mtch_toler)[0]
                # ID
                newID += 1
                IDs[new_idx[isep]] = newID
                dup_filled[isep] = True  # Avoids the other dup(s)
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    '''
    if close_pairs: # TEST SDSS
        tc = SkyCoord(ra=210.053222, dec=31.581701, unit='deg')#, (210.053552, 31.58131)]>
        isep = np.argmin(tc.separation(c_new))
        pdb.set_trace()
        IDs[isep]
    '''
    # Return
    return IDs


def init_data(npix, include_co=False):
    """ Generate an empty masked array for a spectral dataset

    Parameters
    ----------
    npix : int
    include_co : bool, optional
      Include a continuum?

    Returns
    -------
    data : masked  ndarray
    """
    dtype = [(str('wave'), 'float64', (npix)),
          (str('flux'), 'float32', (npix)),
          (str('sig'),  'float32', (npix))]
    if include_co:
        dtype += [(str('co'),   'float32', (npix))]
    data = np.ma.empty((1,), dtype=dtype)
    # Return
    return data


def set_new_ids(maindb, meta, idkey, chk=True, first=False, debug=False, **kwargs):
    """ Set the new IDs

    Parameters
    ----------
    maindb
    meta : Table
    idkey : str
    toler
    chk
    first : bool, optional
      First call to setting the IDs

    Returns
    -------
    cut_db : Table
      Cut to the new sources
    new : bool array
    ids : ID values of newdb
    """
    # IDs
    ids = get_new_ids(maindb, meta, idkey, debug=debug, **kwargs) # Includes new and old
    if debug:
        pdb.set_trace()
    # Crop to rows with new IDs
    if first:
        new = ids >= 0
    else:
        new = ids > 0  # -1 * 0 = 0
    newi = np.where(new)[0]
    # Need unique
    uni, idx_uni = np.unique(ids[newi], return_index=True)
    #
    cut_db = meta[newi[idx_uni]]
    cut_db.add_column(Column(ids[newi[idx_uni]], name=idkey))  # Assumes ordered by ID which is true
    # Reset IDs to all positive
    ids = np.abs(ids)
    meta[idkey] = ids
    #
    return cut_db, new, ids


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
            elif 'COS' in head['INSTRUME']:
                instr = 'COS'
            elif 'ISIS' in head['INSTRUME']:
                instr = 'ISIS'
            elif ('test' in head['INSTRUME']) and ('kp4m' in head['TELESCOP']):  # Kludge for old RCS data
                instr = 'RCS'
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
    elif instr == 'COS':
        try:
            return Rdicts[instr][head['OPT_ELEM'].strip()]
        except KeyError:
            print("Need to add {:s}".format(head['DECKNAME']))
            pdb.set_trace()
    elif instr == 'ISIS':
        try:
            return Rdicts[instr][head['ISIGRAT'].strip()]
        except KeyError:
            print("Need to add {:s}".format(head['ISIGRAT']))
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
            res = Rdicts[instr][head['FILTER3']][head['FPMASK']]
        except KeyError:
            print("Need to add {:s} and/or mask {:s}".format(head['FILTER3'],
                                                             head['FPMASK']))
            pdb.set_trace()
        else:
            return res
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
    elif 'RCS' in instr: # KPNO (retired)
        res_1 = Rdicts[instr][head['DISPERSE']]
        swidth = defs.slit_width(head['APERTURE'])
        return res_1 / swidth
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


def set_sv_idkey(idkey):
    """ Set the idkey in the event that it was not set by
     calling start_maindb

    Parameters
    ----------
    idkey : str

    """
    global sv_idkey
    sv_idkey = idkey

def start_maindb(idkey, **kwargs):
    """ Start the main DB catalog

    Returns
    -------
    maindb : Table
    tkeys : list
      List of columns in the table

    """
    # For checking later
    global sv_idkey
    sv_idkey = idkey
    #
    idict = defs.get_db_table_format(**kwargs)
    tkeys = list(idict.keys())
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)
    # ID_key -- should be unique to the database
    maindb[idkey] = -1  # To get the indexing right
    tkeys += [idkey]
    # Return
    return maindb, tkeys


def write_hdf(hdf, dbname, maindb, zpri, gdict, version, epoch=2000.,
              spaceframe='ICRS', **kwargs):
    """
    Parameters
    ----------
    hdf
    dbname
    maindb
    zpri
    gdict
    version : str
    epoch : float, optional
    spaceframe : str, optional

    Returns
    -------

    """
    import json
    import datetime
    # Write
    clean_table_for_hdf(maindb)
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['NAME'] = str.encode(dbname)
    hdf['catalog'].attrs['EPOCH'] = epoch
    hdf['catalog'].attrs['EQUINOX'] = epoch
    hdf['catalog'].attrs['SpaceFrame'] = str.encode(spaceframe)
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['GROUP_DICT'] = json.dumps(ltu.jsonify(gdict))
    hdf['catalog'].attrs['CREATION_DATE'] = str.encode(datetime.date.today().strftime('%Y-%b-%d'))
    hdf['catalog'].attrs['VERSION'] = str.encode(version)
    # kwargs
    for key in kwargs:
        hdf['catalog'].attrs[str.encode(key)] = kwargs[key]
    # Close
    hdf.close()


