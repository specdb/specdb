""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb

def show_group_meta(meta, meta_keys=None, show_all_keys=True):
    """ Show (nicely) a set of meta data

    Parameters
    ----------
    meta_keys : list, optional
      List of keys to put at the front when showing
    show_all_keys : bool, optional
      Show all keys in meta table

    """
    if meta_keys is None:
        mkeys = []
        # Grab IDs
        for key in meta.keys():
            if '_ID' in key:
                mkeys += [key]
        mkeys += ['RA', 'DEC', 'zem', 'SPEC_FILE']
    else:
        mkeys = meta_keys
    # Add in the rest
    if show_all_keys:
        for key in meta.keys():
            if key not in mkeys:
                mkeys += [key]
    meta[mkeys].pprint(max_width=120)
    return

