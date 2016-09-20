""" Module for DB utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob
import warnings
import pdb

def grab_dbfile():
    """ Grabs the DB file
    Returns
    -------
    db_file : str
      full path to the DB file

    """
    if os.getenv('IGMSPEC_DB') is None:
        warnings.warn('Environmental variable IGMSPEC_DB not set.  Assuming this is a test')
        import igmspec
        db_dir = igmspec.__path__[0]+'/tests/files/'
    else:
        db_dir = os.getenv('IGMSPEC_DB')
    #
    fils = glob.glob(db_dir+'/IGMspec_DB_*hdf5')
    fils.sort()
    db_file = fils[-1]  # Should grab the latest
    # Return
    return db_file
