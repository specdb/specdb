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

    wolverine> specdb_get_uvqs -h
    usage: specdb_get_uvqs [-h] [-v VERSION]

    Grab the UVQS DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate
