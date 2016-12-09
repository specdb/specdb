.. highlight:: igmspec

*******
igmspec
*******

IgmSpec is the first public database available with specdb.
It provides spectra that probe the intergalactic medium
(IGM) from published datasets including SDSS,
2QZ, and data from HST, Keck, VLT, and more.

Versions
========

========  ====================================================
Version   Datasets
========  ====================================================
v01       BOSS_DR12, SDSS_DR7, KODIAQ_DR1, HD_LLS, GGG
v02       MUSoDLA, HSTQSO, COS-Dwarfs, COS-Halos, 2QZ, HDLA100
  ..      ESI-DLA, XQ-100, HST_z2
========  ====================================================

Downloading
===========

Use the script `specdb_get_igmspec` to grab a copy of the database.
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

A greater description can be found in the docs of that repository:
`igmspec <http://github.com/specdb/igmspec>`_.

