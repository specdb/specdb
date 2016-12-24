""" Module for SpecDB Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np
import warnings
import h5py

import datetime
import json

from astropy import units as u

from linetools import utils as ltu

from specdb.specdb import IgmSpec
import specdb



def igmspec_file(version='v02', nspec=5):
    """ Build a debug file from IGMspec
    Returns
    -------

    """
    # Load IGMSpec
    igmsp = IgmSpec(version='v01')  # Will advance as needed
    #
    outfil = specdb.__path__[0]+'/tests/files/IGMspec_DB_{:s}_debug.hdf5'.format(version)
    hdf = h5py.File(outfil,'w')
    # Grab 100 sources from several datasets
    dsets = ['BOSS_DR12', 'HD-LLS_DR1', 'SDSS_DR7', 'GGG']#, '2QZ']
    flags = igmsp.cat['flag_group'].data
    all_IDs = []
    for dset in dsets:
        sflag = igmsp.group_dict[dset]
        query = (flags % (sflag*2)) >= sflag
        gdi = np.where(query)[0]
        # Take nspec
        keep = gdi[0:nspec]
        all_IDs += keep.tolist()  # Save for main catalog
        # Grab data
        if dset == '2QZ':
            pdb.set_trace()
        rows = igmsp[dset].ids_to_allrows(keep)   # Match to all rows
        spec, meta = igmsp[dset].grab_specmeta(rows)  # Grab
        # Group
        grp = hdf.create_group(dset)
        spec_set = hdf[dset].create_dataset('spec', data=spec.data, chunks=True, compression='gzip')
        hdf[dset]['meta'] = meta
    IDs = np.unique(np.array(all_IDs))
    # Catalog
    hdf['catalog'] = igmsp.qcat.cat[IDs]
    hdf['catalog'].attrs['NAME'] = 'igmspec'
    hdf['catalog'].attrs['EPOCH'] = 2000.
    zpri = igmsp.hdf['catalog'].attrs['Z_PRIORITY']
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['VERSION'] = version
    hdf['catalog'].attrs['CREATION_DATE'] = str(datetime.date.today().strftime('%Y-%b-%d'))
    hdfkeys = hdf.keys()
    sdict = igmsp.group_dict.copy()
    for dkey in sdict.keys():
        if dkey not in hdfkeys:
            sdict.pop(dkey, None)
    hdf['catalog'].attrs['GROUP_DICT'] = json.dumps(ltu.jsonify(sdict))
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))

def main(flg_file, sdss=None, ml_survey=None):

    # Sightlines
    if (flg_file % 2**1) >= 2**0:
        igmspec_file()



# Run
if __name__ == '__main__':
    # Run from above src/
    #  I.e.   python src/training_set.py
    flg_file = 0
    flg_file += 2**0   # IGMspec

    main(flg_file)
