""" Module for DB utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob
import warnings
import pdb


def load_db(db_type, **kwargs):
    """
    Parameters
    ----------
    db_type : str
      Type of database
      Current allowed entries are [igmspec]


    Returns
    -------
    dbobj : SpecDB object
      Generally a child

    """
    if db_type == 'igmspec':
        from specdb.specdb import IgmSpec
        Specdb = IgmSpec(**kwargs)
    else:
        raise IOError("Not ready for this dbase value: {:s}".format(args.dbase))

    # Return
    return Specdb

