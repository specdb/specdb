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
v06       All of the spectra of the QPQI-X series
v07       Coming soon;  Includes near-IR spectra from Coatman+18
========  ======================================================

Catalog
=======

The source catalog includes the QPQX catalog published by
Finlay et al. 2018 as the first 3836 entries.  Additional
meta data is provided for these sources (e.g. GALEX, WISE fluxes).
A significant number, however, have no spectra within this
specDB database.  Most are from BOSS/SDSS and may
be found in the igmspec specDB.  The remainder were discovered
at APO by JFH and are not available.

There are also approximately 1800 spectra included in the specDB
that are not listed in the QPQX catalog of Findlay et al. 2018.
Nearly all of these are BOSS and SDSS spectra that were continuum
fit as part of the QPQ survey and we are archiving these continua.

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
