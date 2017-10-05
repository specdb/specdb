.. highlight:: rest

**************
SDSS in SpecDB
**************

Because SDSS/BOSS is likely to be a major contributor
to many specDB files, there are a few special methods
and scripts for those data.


get_sdss
========

A default method of the `specdb-class`_ is
get_sdss() which takes plate, fiber and returns
meta and spectra::

    spec, meta = Specdb.get_sdss(plate, fiberid)

specdb_sdss
===========

This script also takes plate, fiber and shows
meta data and may plot the spectra.  See
`sdss-spec`_ for more details.

