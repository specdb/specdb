""" Module to build a private DB
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import json
import h5py
import warnings
import pdb

from astropy.table import Table, vstack, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.time import Time

from linetools import utils as ltu
from linetools.spectra import io as lsio
from linetools.spectra.xspectrum1d import XSpectrum1D

from specdb.build import utils as spbu
from specdb.zem import utils as spzu
from specdb import defs


def grab_files(tree_root, skip_files=('c.fits', 'C.fits', 'e.fits',
                                      'E.fits', 'N.fits', 'old.fits'),
               only_conti=False):
    """ Generate a list of FITS files within the file tree

    Parameters
    ----------
    tree_root : str
      Top level path of the tree of FITS files
    skip_files : tuple
      List of file roots to skip as primary files when ingesting
    only_conti : bool, optional
      Only grab files with separate continua files (mainly for QPQ)

    Returns
    -------
    files : list
      List of FITS files

    """
    walk = os.walk(tree_root)
    folders = ['/.']
    pfiles = []
    while len(folders) > 0:
        # Search for fits files
        ofiles = []
        for folder in folders:
            if only_conti:
                ofiles += glob.glob(tree_root+'/'+folder+'/*_c.fits*')
            else:
                ofiles += glob.glob(tree_root+'/'+folder+'/*.fits*')
            # Eliminate error and continua files
            for ofile in ofiles:
                flg = True
                # Ugly loop
                for skip_file in skip_files:
                    if skip_file in ofile:
                        flg = False
                #
                if only_conti:
                    ofile = ofile.replace('_c','')
                    if not os.path.isfile(ofile):
                        print("{:s} not present".format(ofile))
                if flg:
                    pfiles.append(ofile)
            # walk
        folders = next(walk)[1]
    # Return
    return pfiles


def mk_meta(files, ztbl, fname=False, stype='QSO', skip_badz=False,
            mdict=None, parse_head=None, debug=False, chkz=False,
            verbose=False, specdb=None, sdb_key=None, **kwargs):
    """ Generate a meta Table from an input list of files

    Parameters
    ----------
    files : list
      List of FITS files
    ztbl : Table
      Table of redshifts.  Must include RA, DEC, ZEM, ZEM_SOURCE
      Used for RA/DEC if fname=False;  then requires SPEC_FILE too
    fname : bool, optional
      Attempt to parse RA/DEC from the file name
      Format must be
      SDSSJ######(.##)+/-######(.#)[x]
        where x cannot be a #. or +/-
    specdb : SpecDB, optional
      Database object to grab ID values from
      Requires sdb_key
    sdb_key : str, optional
      ID key in SpecDB object
    skip_badz : bool, optional
      Skip spectra without a parseable redshift (using the Myers catalog)
    parse_head : dict, optional
      Parse header for meta info with this dict
    mdict : dict, optional
      Input meta data in dict form e.g.  mdict=dict(INSTR='ESI')
    chkz : bool, optional
      If any sources have no parseable redshift, hit a set_trace

    Returns
    -------
    meta : Table
      Meta table
    """
    if specdb is not None:
        if sdb_key is None:
            raise IOError("Must specify sdb_key if you are passing in specdb")
    Rdicts = defs.get_res_dicts()
    #
    coordlist = []
    for ifile in files:
        if fname:
            # Starting index
            if 'SDSSJ' in ifile:
                i0 = ifile.find('SDSSJ')+4
            else:
                i0 = ifile.rfind('J')+1
            # Find end (ugly)
            for ii in range(i0+1,99999):
                if ifile[ii] in ('0','1','2','3','4','5','6','7','8','9',
                                 '.','+','-'):
                    continue
                else:
                    i1 = ii
                    break
            # Deal with .fits
            if ifile[i1-1] == '.':
                i1 -= 1
            # Get coord
            try:
                coord = ltu.radec_to_coord(ifile[i0:i1])
            except (UnboundLocalError, ValueError):
                pdb.set_trace()
        else:
            sname = ifile.split('/')[-1]
            mt = np.where(ztbl['SPEC_FILE'] == sname)[0]
            if len(mt) != 1:
                raise IndexError("NO MATCH FOR {:s}".format(sname))
            coord = ltu.radec_to_coord((ztbl['RA'][mt],
                                        ztbl['DEC'][mt]))[0]
        coordlist.append(coord)
    coords = SkyCoord(ra=[coord.ra.degree for coord in coordlist], dec=[coord.dec.degree for coord in coordlist], unit='deg')

    # Generate Meta Table
    maindb, tkeys = spbu.start_maindb(private=True)

    # Fill
    meta = Table()
    meta['RA'] = coords.ra.deg
    meta['DEC'] = coords.dec.deg
    meta['STYPE'] = [stype]*len(meta)
    meta['flag_survey'] = [1]*len(meta)

    zem, zsource = spzu.zem_from_radec(meta['RA'], meta['DEC'], ztbl)
    badz = zem <= 0.
    if np.sum(badz) > 0:
        if skip_badz:
            warnings.warn("Skipping {:d} entries without a parseable redshift".format(
                np.sum(badz)))
        else:
            # Try specdb if it exists
            if specdb is not None:
                pdb.set_trace()
                sztbl = specdb.qcat.cat[['RA','zem','DEC', 'flag_zem']]
                sztbl.rename_column('flag_zem', 'ZEM_SOURCE')
                sztbl.rename_column('zem', 'ZEM')
                szem, szsource = spzu.zem_from_radec(meta['RA'][badz], meta['DEC'][badz], sztbl)
                pdb.set_trace()
            if chkz:  # Turn this on to hit a stop instead of an Exception
                pdb.set_trace()
            else:
                raise ValueError("{:d} entries without a parseable redshift".format(
                    np.sum(badz)))
    meta['zem'] = zem
    pdb.set_trace()
    meta['sig_zem'] = 0.  # Need to add
    meta['flag_zem'] = zsource
    # Cut
    meta = meta[~badz]

    # specdb IDs
    if specdb is not None:
        meta[sdb_key] = [-9999]*len(meta)
        c_igmsp = SkyCoord(ra=specdb.qcat.cat['RA'], dec=specdb.qcat.cat['DEC'], unit='deg')
        c_new = SkyCoord(ra=meta['RA'], dec=meta['DEC'], unit='deg')
        # Find new sources
        idx, d2d, d3d = match_coordinates_sky(c_new, c_igmsp, nthneighbor=1)
        cdict = defs.get_cat_dict()
        mtch = d2d < cdict['match_toler']
        meta[sdb_key][mtch] = specdb.qcat.cat[sdb_key][idx[mtch]]

    # Stack (primarily as a test)
    try:
        maindb = vstack([maindb,meta], join_type='exact')
    except:
        pdb.set_trace()
    maindb = maindb[1:]

    # SPEC_FILE
    maindb['SPEC_FILE'] = np.array(files)[~badz]

    # Try Header?
    if parse_head is not None:
        # Setup to store
        plist = {}
        for key in parse_head.keys():
            plist[key] = []
        # Loop on files
        for sfile in maindb['SPEC_FILE']:
            if verbose:
                print('Parsing {:s}'.format(sfile))
            head = fits.open(sfile)[0].header
            for key,item in parse_head.items():
                # R
                if key == 'R':
                    if parse_head[key] == True:
                        try:
                            plist[key].append(spbu.set_resolution(head))
                        except ValueError:
                            if mdict is not None:
                                try:
                                    plist[key].append(mdict['R'])
                                except KeyError:
                                    pdb.set_trace()
                            else:
                                pdb.set_trace()
                                plist[key].append(0.)
                    else:
                        raise ValueError("Set something else for R")
                elif key == 'DATE-OBS':
                    tval = Time(head[item].replace('/','-'), format='isot', out_subfmt='date')
                    plist[key].append(tval.iso)
                else:
                    plist[key].append(head[item])
            # INSTRUMENT SPECIFIC
            try:
                instr = head['INSTRUME']
            except KeyError:
                instr = 'none'
            if 'LRIS' in instr:
                if 'GRATING' not in plist.keys():
                    plist['GRATING'] = []
                    plist['INSTR'] = []
                    plist['R'] = []
                try:
                    det = head['DETECTOR']
                except KeyError:
                    if head['OUTFILE'] == 'lred':
                        det = 'LRIS-R'
                    else:
                        det = 'LRIS-B'
                if 'LRIS-R' in det:
                    plist['GRATING'].append(head['GRANAME'])
                    plist['INSTR'].append('LRISr')
                else:
                    plist['GRATING'].append(head['GRISNAME'])
                    plist['INSTR'].append('LRISb')
                # Resolution
                res = Rdicts[plist['INSTR'][-1]][plist['GRATING'][-1]]
                try:
                    sname = head['SLITNAME']
                except KeyError:
                    swidth = 1.
                else:
                    swidth = defs.slit_width(sname)
                plist['R'].append(res/swidth)
        # Finish
        for key in plist.keys():
            maindb[key] = plist[key]
    # mdict
    if mdict is not None:
        for key,item in mdict.items():
            maindb[key] = [item]*len(meta)

    # EPOCH
    if 'EPOCH' not in maindb.keys():
        warnings.warn("EPOCH not defined.  Filling with 2000.")
        maindb['EPOCH'] = 2000.

    # Fill in empty columns with warning
    mkeys = maindb.keys()
    req_clms = defs.get_req_clms()
    for clm in req_clms:
        if clm not in mkeys:
            if clm not in ['NPIX','WV_MIN','WV_MAX']:  # File in ingest_spec
                warnings.warn("Meta Column {:s} not defined.  Filling with DUMMY".format(clm))
                if clm == 'DATE-OBS':
                    maindb[clm] = ['9999-1-1']*len(maindb)
                else:
                    maindb[clm] = ['DUMMY']*len(maindb)

    # Return
    if debug:
        maindb[['IGM_ID', 'RA', 'DEC', 'SPEC_FILE']].pprint(max_width=120)
        pdb.set_trace()
    return maindb


def dumb_spec():
    """ Generate a dummy spectrum
    Returns
    -------

    """
    npix = 1000
    dspec = XSpectrum1D.from_tuple((np.arange(npix)+5000., np.ones(npix),
                                   np.ones(npix)))
    #
    return dspec


def ingest_spectra(hdf, sname, meta, max_npix=10000, chk_meta_only=False,
                   refs=None, verbose=False, badf=None, grab_conti=False, **kwargs):
    """ Ingest the spectra
    Parameters
    ----------
    hdf : hdf5 pointer
    sname : str
      Name of dataset
    meta : Table
    max_npix : int, optional
      Maximum length of the spectra
    chk_meta_only : bool, optional
      Only check meta file;  will not write
    refs : list, optional
      list of dicts with reference info
    badf : list, optional
      List of bad spectra [use only if you know what you are doing!]
    grab_conti : bool, optional
      Grab continua.  They should exist but do not have to

    Returns
    -------

    """
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    grp = hdf.create_group(sname)
    # Spectra
    nspec = len(meta)
    dtypes=[(str('wave'), 'float64', (max_npix)),
           (str('flux'), 'float32', (max_npix)),
           (str('sig'),  'float32', (max_npix))]
    dkeys = ['wave','flux','sig']
    if grab_conti:
        dtypes += [(str('co'),   'float32', (max_npix))]
        dkeys += ['co']
    data = np.ma.empty((1,), dtype=dtypes)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    # Loop
    for jj,member in enumerate(meta['SPEC_FILE']):
        # Extract
        f = member
        # Parse name
        fname = f.split('/')[-1]
        if verbose:
            print(fname)
        # Read
        if badf is not None:
            for ibadf in badf:
                if ibadf in f:
                    spec = dumb_spec()
                else:
                    spec = lsio.readspec(f)
        else:
            spec = lsio.readspec(f)
        # npix
        head = spec.header
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        # Some fiddling about
        for key in dkeys:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        if grab_conti:
            if spec.co_is_set:
                data['co'][0][:npix] = spec.co.value
        # Meta
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        # Set
        spec_set[jj] = data

    # Add columns
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))

    # Add HDLLS meta to hdf5
    if spbu.chk_meta(meta):#, skip_igmid=True):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        pdb.set_trace()
        raise ValueError("meta file failed")
    # References
    if refs is not None:
        jrefs = ltu.jsonify(refs)
        hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


def mk_db(trees, names, outfil, ztbl, **kwargs):
    """ Generate the DB

    Parameters
    ----------
    trees : list
      List of top level paths for the FITS files
    names : list
      List of names for the various datasets
    outfil : str
      Output file name for the hdf5 file
    ztbl : Table
      See above

    Returns
    -------

    """
    from igmspec import defs as igmsp_defs
    # HDF5 file
    hdf = h5py.File(outfil,'w')

    # Defs
    zpri = igmsp_defs.z_priority()
    sdict = {}

    # Main DB Table
    maindb, tkeys = spbu.start_maindb(private=True)
    maindb['PRIV_ID'] = -1  # To get the indexing right
    tkeys += ['PRIV_ID']

    # MAIN LOOP
    for ss,tree in enumerate(trees):
        print('Working on tree: {:s}'.format(tree))
        # Files
        fits_files = grab_files(tree)
        # Meta
        full_meta = mk_meta(fits_files, ztbl, **kwargs)
        # Survey IDs
        flag_s = 2**ss
        sdict[names[ss]] = flag_s
        if ss == 0:
            ids = np.arange(len(full_meta), dtype=int)
            full_meta['PRIV_ID'] = ids
            full_meta['flag_survey'] = flag_s
            cut = full_meta
        else:
            cut, new, ids = spbu.set_new_ids(maindb, full_meta, idkey='PRIV_ID')
            cut['flag_survey'] = [flag_s]*len(cut)
            midx = np.array(maindb['PRIV_ID'][ids[~new]])
            maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
        # Catalog
        cat_meta = cut[tkeys]
        assert spbu.chk_maindb_join(maindb, cat_meta)
        # Append
        maindb = vstack([maindb,cat_meta], join_type='exact')
        if ss == 0:
            maindb = maindb[1:]  # Eliminate dummy line
        # Ingest
        ingest_spectra(hdf, names[ss], full_meta, **kwargs)

    # Write
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['EPOCH'] = 2000.
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['SURVEY_DICT'] = json.dumps(ltu.jsonify(sdict))
    #hdf['catalog'].attrs['VERSION'] = version
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    #hdf['catalog'].attrs['SURVEY_DICT'] = defs.get_survey_dict()
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))


