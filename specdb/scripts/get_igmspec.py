#!/usr/bin/env python
"""
Run a build of the DB
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Grab the IGMspec DB')
    parser.add_argument("-v", "--version", default='v02', help="DB version to generate")
    #parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs):
    """ Run
    Parameters
    ----------
    pargs

    Returns
    -------

    """
    import subprocess

    # Version
    if pargs.version not in ['v01', 'v02']:
        raise IOError("Bad version number")

    # URL
    #url_page = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/'
    url_page = 'http://specdb.ucsc.edu/'
    url = url_page+'IGMspec_DB_{:s}.hdf5'.format(pargs.version)
    
    # wget command
    subprocess.call(['wget', '--continue', '--timestamping', url])

if __name__ == '__main__':
    # Check for wget
    if not any(os.access(os.path.join(path, 'wget'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
        raise RuntimeError("You need to install wget in your PATH")
    # Giddy up
    main()
