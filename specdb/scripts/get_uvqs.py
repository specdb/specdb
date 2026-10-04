#!/usr/bin/env python
"""
Grab the UVQS DB
"""

import os


def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Grab the UVQS DB')
    parser.add_argument("-v", "--version", default='v01', help="DB version to generate")
    #parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(pargs=None):
    """ Run
    Parameters
    ----------
    pargs

    Returns
    -------

    """
    if pargs is None:
        pargs = parser()
    import subprocess

    # Version
    if pargs.version not in ['v01', 'v02']:
        raise IOError("Bad version number")

    # URL
    url_page = 'http://specdb.ucsc.edu/'
    url = url_page+'UVQS_DB_{:s}.hdf5'.format(pargs.version)
    
    # wget command
    subprocess.call(['wget', '--continue', '--timestamping', url])

if __name__ == '__main__':
    # Check for wget
    if not any(os.access(os.path.join(path, 'wget'), os.X_OK) for path in os.environ["PATH"].split(os.pathsep)):
        raise RuntimeError("You need to install wget in your PATH")
    # Giddy up
    main()
