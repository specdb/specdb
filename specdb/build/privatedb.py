""" Module to build a private DB
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import json
import h5py
import warnings
import pdb
import datetime

from astropy.table import Table, Column
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.time import Time

from linetools import utils as ltu
from linetools.spectra import io as lsio
from linetools.spectra.xspectrum1d import XSpectrum1D

from specdb.build import utils as spbu
from specdb.zem import utils as spzu
from specdb import defs
from specdb.build.utils import add_ids, write_hdf

try:
    basestring
except NameError:  # For Python 3
    basestring = str


def grab_files(branch, skip_files=('c.fits', 'C.fits', 'e.fits',
                                      'E.fits', 'N.fits', 'old.fits'),
               only_conti=False, skip_folders=[]):
    """ Generate a list of FITS files within the file tree

    Parameters
    ----------
    branch : str
      Top level path of the FITS files
    skip_files : tuple
      List of file roots to skip as primary files when ingesting
    only_conti : bool, optional
      Only grab files with separate continua files (mainly for QPQ)
    skip_folders : list, optional
      Skip any folder with these names

    Returns
    -------
    pfiles : list
      List of FITS files
    meta_file : str or None
      Name of meta file in tree_root
    """
    walk = os.walk(branch)
    folders = ['/.']
    pfiles = []
    while len(folders) > 0:
        # Search for fits files
        ofiles = []
        for folder in folders:
            if folder in skip_folders:
                print("Skipping folder = {:s}".format(folder))
                continue
            if only_conti:
                ofiles += glob.glob(branch+'/'+folder+'/*_c.fits*')
            else:
                ofiles += glob.glob(branch+'/'+folder+'/*.fits*')
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
    # Grab meta file (if one exists)
    mfile = glob.glob(branch+'/*_meta.json')
    if len(mfile) == 1:
        mfile = mfile[0]
    elif len(mfile) == 0:
        warnings.warn("No meta JSON file found in branch: {:s}  This is not recommended".format(branch))
        mfile = None
    else:
        raise IOError("Multiple meta JSON files in branch: {:s}.  Limit to one".format(branch))
    # Return
    return pfiles, mfile


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
    stype : str, optional
      Description of object type (e.g. 'QSO', 'Galaxy', 'SN')
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
    ras = np.array([coord.ra.degree for coord in coordlist])
    decs= np.array([coord.dec.degree for coord in coordlist])
    coords = SkyCoord(ra=ras, dec=decs, unit='deg')

    # Generate maindb Table
    #maindb, tkeys = spbu.start_maindb(private=True)

    # Fill
    meta = Table()
    meta['RA_GROUP'] = coords.ra.deg
    meta['DEC_GROUP'] = coords.dec.deg
    meta['STYPE'] = [str(stype)]*len(meta)

    zem, zsource = spzu.zem_from_radec(meta['RA_GROUP'], meta['DEC_GROUP'], ztbl)
    badz = zem <= 0.
    if np.sum(badz) > 0:
        if skip_badz:
            warnings.warn("Skipping {:d} entries without a parseable redshift".format(
                np.sum(badz)))
        else:
            if chkz:  # Turn this on to hit a stop instead of an Exception
                pdb.set_trace()
            else:
                raise ValueError("{:d} entries without a parseable redshift".format(
                    np.sum(badz)))
    meta['zem'] = zem
    meta['sig_zem'] = 0.  # Need to add
    meta['flag_zem'] = zsource
    # Cut
    meta = meta[~badz]

    # specdb IDs
    if sdb_key is not None:
        meta[sdb_key] = [-9999]*len(meta)
        if sdb_key not in meta.keys():
            meta[sdb_key] = [-9999]*len(meta)
        c_igmsp = SkyCoord(ra=specdb.qcat.cat['RA'], dec=specdb.qcat.cat['DEC'], unit='deg')
        c_new = SkyCoord(ra=meta['RA_GROUP'], dec=meta['DEC_GROUP'], unit='deg')
        # Find new sources
        idx, d2d, d3d = match_coordinates_sky(c_new, c_igmsp, nthneighbor=1)
        cdict = defs.get_cat_dict()
        mtch = d2d < cdict['match_toler']
        meta[sdb_key][mtch] = specdb.qcat.cat[sdb_key][idx[mtch]]

    # Stack (primarily as a test)
    '''
    try:
        maindb = vstack([maindb,meta], join_type='exact')
    except:
        pdb.set_trace()
    '''

    # SPEC_FILE
    meta['SPEC_FILE'] = np.array(files)[~badz]

    # Try Header?
    if parse_head is not None:
        # Setup to store
        plist = {}
        for key in parse_head.keys():
            plist[key] = []
        # Loop on files
        for sfile in meta['SPEC_FILE']:
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
            meta[key] = plist[key]
    # mdict
    if mdict is not None:
        for key,item in mdict.items():
            meta[key] = [item]*len(meta)

    # EPOCH
    if 'EPOCH' not in meta.keys():
        warnings.warn("EPOCH not defined.  Filling with 2000.")
        meta['EPOCH'] = 2000.

    # GROUP ID
    meta['GROUP_ID'] = np.arange(len(meta)).astype(int)

    # Fill in empty columns with warning
    mkeys = meta.keys()
    req_clms = defs.get_req_clms(sdb_key=sdb_key)
    for clm in req_clms:
        if clm not in mkeys:
            if clm not in ['NPIX','WV_MIN','WV_MAX']:  # File in ingest_spec
                warnings.warn("Meta Column {:s} not defined.  Filling with DUMMY".format(clm))
                if clm == 'DATE-OBS':
                    meta[clm] = ['9999-1-1']*len(meta)
                else:
                    meta[clm] = ['DUMMY']*len(meta)

    # Return
    if debug:
        meta[['RA_GROUP', 'DEC_GROUP', 'SPEC_FILE']].pprint(max_width=120)
        pdb.set_trace()
    return meta


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
                   refs=None, verbose=False, badf=None,
                   grab_conti=False, **kwargs):
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
    print("Adding {:s} group to DB".format(sname))
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
                    try:
                        spec = lsio.readspec(f)#, **kwargs)
                    except ValueError:  # Probably a continuum problem
                        pdb.set_trace()
        else:
            spec = lsio.readspec(f)#, **kwargs)
        # npix
        head = spec.header
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d} vs {:d})".format(
                    npix, max_npix))
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


def mk_db(dbname, tree, outfil, iztbl, version='v00', **kwargs):
    """ Generate the DB

    Parameters
    ----------
    dbname : str
      Name for the database
    tree : str
      Path to top level of the tree of FITS files
      Typically, each branch in the tree corresponds to a single instrument
    outfil : str
      Output file name for the hdf5 file
    iztbl : Table or str
      If Table, see meta() docs for details on its format
      If str, it must be 'igmspec' and the user must have that DB downloaded
    version : str, optional
      Version code

    Returns
    -------

    """
    from specdb import defs

    # ztbl
    if isinstance(iztbl, basestring):
        if iztbl == 'igmspec':
            from specdb.specdb import IgmSpec
            igmsp = IgmSpec()
            ztbl = Table(igmsp.idb.hdf['quasars'].value)
    elif isinstance(iztbl, Table):
        ztbl = iztbl
    else:
        raise IOError("Bad type for ztbl")

    # Find the branches
    branches = glob.glob(tree+'/*')
    branches.sort()
    # HDF5 file
    hdf = h5py.File(outfil,'w')

    # Defs
    zpri = defs.z_priority()
    gdict = {}

    # Main DB Table
    maindb, tkeys = spbu.start_maindb('PRIV_ID')

    # MAIN LOOP
    for ss,branch in enumerate(branches):
        # Skip files
        if not os.path.isdir(branch):
            continue
        print('Working on branch: {:s}'.format(branch))
        # Files
        fits_files, meta_file = grab_files(branch)
        # Meta
        maxpix, phead, mdict, stype = 10000, None, None, 'QSO'
        if meta_file is not None:
            # Load
            meta_dict = ltu.loadjson(meta_file)
            # Maxpix
            if 'maxpix' in meta_dict.keys():
                maxpix = meta_dict['maxpix']
            # STYPE
            if 'stype' in meta_dict.keys():
                stype = meta_dict['stype']
            # Parse header
            if 'parse_head' in meta_dict.keys():
                phead = meta_dict['parse_head']
            if 'meta_dict' in meta_dict.keys():
                mdict = meta_dict['meta_dict']
        full_meta = mk_meta(fits_files, ztbl,
                            parse_head=phead, mdict=mdict, **kwargs)
        # Update group dict
        group_name = branch.split('/')[-1]
        flag_g = spbu.add_to_group_dict(group_name, gdict)
        # IDs
        maindb = add_ids(maindb, full_meta, flag_g, tkeys, 'PRIV_ID', first=(flag_g==1))
        # Ingest
        ingest_spectra(hdf, group_name, full_meta, max_npix=maxpix, **kwargs)

    # Write
    write_hdf(hdf, dbname, maindb, zpri, gdict, version)
    print("Wrote {:s} DB file".format(outfil))


