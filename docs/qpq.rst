.. highlight:: qpq

***
qpq
***

QPQ is the third public database available with specdb.
It provides all of the spectra from the QPQ survey
(Findlay et al. 2018).

Versions
========

The dataset(s) included in a specific version of the
`qpq` database release are summarized below.

========  ======================================================
Version   Data Groups
========  ======================================================
v05       All of the spectra of the QPQI-X series
v06       Coming soon;  Includes near-IR spectra from Coatman+18
========  ======================================================

Downloading
===========

Use the script `specdb_get_qpq` to grab a copy of the database.
Here is the usage::

    profx> specdb_get_qpq -h
    usage: specdb_get_qpq [-h] [-v VERSION]

    Grab the QPQ DB [v1.0]

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to grab [v05]

CAT_META
========

In addition to the main catalog which provides brief all the sources
