""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings

def show_group_meta(meta, meta_keys=None, show_all_keys=True, idkey=None, **kwargs):
    """ Show (nicely) a set of meta data

    Parameters
    ----------
    meta_keys : list, optional
      List of keys to put at the front when showing
    show_all_keys : bool, optional
      Show all keys in meta table
    **kwargs : Passed to pprint() of the Table

    """
    if meta_keys is None:
        mkeys = []
        # Grab IDs
        if idkey is None:
            for key in meta.keys():
                if '_ID' in key:
                    mkeys += [key]
        else:
            mkeys += [idkey]
        mkeys += ['GROUP', 'RA_GROUP', 'DEC_GROUP', 'zem_GROUP', 'SPEC_FILE']
    else:
        mkeys = meta_keys
    # Add in the rest
    if show_all_keys:
        for key in meta.keys():
            if key not in mkeys:
                mkeys += [key]
    # Confirm the keys are in meta
    keep_keys = []
    for key in mkeys:
        if key in meta.keys():
            keep_keys += [key]
        else:
            warnings.warn("Key: {:s} not in meta so not showing".format(key))
    meta[keep_keys].pprint(max_width=120, **kwargs)
    return

