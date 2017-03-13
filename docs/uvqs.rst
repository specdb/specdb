.. highlight:: uvqs

****
uvqs
****

UVQS is the second public database available with specdb.
It provides all of the spectra from the UVQS survey.  These
are also available at MAST.

Versions
========

The dataset(s) included in a specific version of the
`uvqs` database release are summarized below.

========  ====================================================
Version   Data Groups
========  ====================================================
v01       FUV
v02       To be released with NUV sample
========  ====================================================

Downloading
===========

Use the script `specdb_get_uvqs` to grab a copy of the database.
Here is the usage::

    usage: specdb_get_igmspec [-h] [-v VERSION]

    Grab the IGMspec DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate


The standard call is therefore::

    UNIX> cd where_you_want_the_DB_file
    UNIX> specdb_get_igmspec -v=v02

Setup
=====

Add to your UNIX environmental variables the path to the DB file, e.g.::

    setenv IGMSPEC_DB /raid/IGMSPEC_DB/   # csh or tcsh
    export IGMSPEC_DB="/raid/IGMSPEC_DB/" # bash


More
====

You can also view a series of examples of using the database here:
`Examples with igmspec <https://github.com/specdb/specdb/blob/master/docs/nb/Examples_with_igmspec.ipynb>`_

Other details are provided in the
`igmspec documentation <http://igmspec.readthedocs.io/en/latest/>`_.

