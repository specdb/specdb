.. highlight:: catalog

*******
Catalog
*******

This document describes the source catalog of
a `specdb` database
and methods to query it.

.. _catalog-desc:

Description
===========

Fundamental to a `specdb` database is the source catalog, written
to the 'catalog' dataset in the hdf5 file.  The catalog table
contains the complete list of sources with at least one spectrum
in the `specdb` database.  It is intended that each entry refers
to a unique source on the sky.  Typically, each source is a unique
physical object although there are obvious exceptions, e.g. multiple
lens images of a quasar.

At a minimum the catalog
table contains the RA, DEC, redshift info, a group flag indicating
which data groups contain spectra on the source,
and a unique catalog IDKEY identifier.

The table comprising the source catalog has the following entries:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IDKEY       int      Unique identifier;  name of the key should include ID (e.g. IGM_ID)
flag_group  int      Bitwise flag indicating the groups that the source has spectra in
zem         float    Emission redshift of background source
sig_zem     float    Estimated error in the redshift
flag_zem    str      Key indicating source of the redshift (e.g. BOSS_PCA)
RA          float    Right Ascension (deg)
DEC         float    Declination (deg)
STYPE       str      Type of the source (e.g. QSO)
==========  ======== ============================================


.. _query-catalog:

Querying the Source Catalog
===========================

There are several approaches to querying the source catalog
which differ according to the treatment of coordinates and
source IDs.  The examples in the following documentation
work on the test database file provided with `specdb`.
To follow along, instantiate a SpecDB class::

    import specdb
    db_file = specdb.__path__[0]+'/tests/files/IGMspec_DB_v02_debug.hdf5'
    from specdb.specdb import SpecDB
    sdb = SpecDB(db_file=db_file)

One can proceed to querying.

Querying without Coordinates
----------------------------

One can query the source catalog without regard
to coordinates.  In this case, a :doc:`query_dict`
is required.  Here is an example::

    qdict = {'zem': (3.,5.), 'flag_group-BITWISE-OR': [2,4,32], 'STYPE': 'QSO'}
    matches, sub_cat, IDs = sdb.qcat.query_dict(qdict)

This query restricts to sources with redshifts 3<z<5,
source type 'QSO', and with spectra in any of the data
groups specified by the bitwise flags of 2,4, or 32.

The method returns a bool array (matches) indicating which
rows of the catalog are matched, the sub-catalog of those rows,
and the IDKEY values of those rows.

Querying a Position on the Sky
------------------------------

Querying with a List of Coordinates
-----------------------------------

There are several methods that interface with the primary
source catalog.

cat_from_coords
---------------

This method returns a Table drawn from the catalog matching
the size and order of an input set of coordinates.  Sources
that are not matched within the tolerance (default = 0.5 arcsec)
have entries filled with zero values and ID<0.

Here is an example call::

    coords = SkyCoord(ra=[0.0028,0.0019], dec=[14.9747,17.7737], unit='deg')
    sub_cat = igmsp.qcat.cat_from_coords(coords)

The user can then analyze the catalog for this subset of
sources (if any matched).

radial_search
-------------

One may search to within a given radius for sources around
an input coordinate.  Here is an example::

   ids2334 = igmsp.qcat.radial_search('J233446.40-090812.3', 1.*u.arcsec)

The method returns an array of all source IDs within that radius.

match_coord
-----------

This method matches a set of input coordinates (a SkyCoord object)
to the source catalog within an optional tolerance (default=0.5").  It returns
an ndarray of IDs with shape and order matching the input list.
Coordinates without a match within the tolerance
have -1 values .  Here is an example::

    coords = SkyCoord(ra=[0.0019,1.2321], dec=[17.7737,-12.2332], unit='deg')
    IDs = igmsp.qcat.match_coord(coords)

One can further restrict the search to a specific group

    IDs = igmsp.match_coord(coords, group='BOSS_DR12')

Sources that are a match in position but not within the group
have an ID=-2.

I/O
===

show
----

A printout of the catalog values for a list of IDs is provided
by `show_cat`::

   igmsp.qcat.show_cat(IDs)

This includes the flag_group values which indicate the groups
that include a given source.  The catalog only shows a single
entry per source and only those sources with ID values within
the catalog (e.g. negative values are ignored).
