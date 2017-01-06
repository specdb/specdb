.. highlight:: meta

*********
Meta data
*********

This document describes the meta data of
a `specdb` database and methods to query it.

.. _meta-desc:

Description
===========

Each of the :doc:`groups` in a `specdb` database
contains a meta data table.


At a minimum the meta table contains the following
columns:

The table comprising the source catalog has the following entries:

    req_clms = ['RA_GROUP', 'DEC_GROUP', 'EPOCH', 'zem_GROUP', 'R', 'WV_MIN',
            'WV_MAX', 'DATE-OBS', 'GROUP_ID', 'NPIX', 'SPEC_FILE',
            'INSTR', 'DISPERSER', 'TELESCOPE']

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IDKEY       int      Unique identifier in the `specdb` database
zem_GROUP   float    Emission redshift of source given by the dataset
RA_GROUP    float    Right Ascension (deg) given by the dataset
DEC_GROUP   float    Declination (deg) given by the dataset
EPOCH       float    Year of epoch
R           float    Spectral resolution (dlambda/lambda; FWHM)
WV_MIN      float    Minimum wavelength value
WV_MAX      float    Maximum wavelength value
NPIX        int      Number of pixels in the spectrum
SPEC_FILE   str      Individual filename of the spectrum
INSTR       str      Instrument used: see specdb.defs.instruments
DISPERSER   str      Dispersing element
TELESCOPE   str      Name of the telescope
==========  ======== ============================================


.. _access-meta:

Accessing the Meta Data
=======================

Within the `specdb` software, a meta data
table for a given group is read into memory
during instantiation of an :ref:`interface-group`
object.  The meta data is stored as
an astropy.table.Table.

For convenience, InterfaceGroup objects
are kept in a hidden dict within the
SpecDB class.

Here is explanation by example::

    # Instantiate the SpecDB object
    sdb = SpecDB(db_file='filename_of_DB')
    # Instantiate the InterfaceGroup and pass a pointer to the meta table
    meta_tbl = sdb['group_name'].meta
    # Accessing data in the meta table
    instruments = meta_tbl['INSTR']

.. _query-meta:

Querying the Meta Data
======================

There are several approaches to querying the meta data
with the SpecDB object.  Each of
these differ according to the treatment of coordinates.
And each uses the underlying :ref:`group-query-meta` method
in the :ref:`interface-group` class.

The examples in the following documentation
work on the test database file provided with `specdb`.
To follow along, instantiate a SpecDB class::

    import specdb
    db_file = specdb.__path__[0]+'/tests/files/IGMspec_DB_v02_debug.hdf5'
    from specdb.specdb import SpecDB
    sdb = SpecDB(db_file=db_file)

One can proceed to querying.
Much of the following is also contained in this
`Query Meta Data Notebook <https://github.com/specdb/specdb/blob/master/docs/nb/Query_Meta.ipynb>`_.


Meta Data with a Query dict
---------------------------

One can query meta data tables with
a :doc:`query_dict` object.
Here is an example::

    qdict = {'TELESCOPE': 'Gemini-North', 'NPIX': (1580,1583), 'DISPERSER': ['B600', 'R400']}
    qmeta = sdb.query_meta(qdict)


This query requests the meta data of all
spectra taken at the Gemini-North telescope
with the GMOS-N instrument
having between 1580-1583 pixels (inclusive)
and with either the B600 or R400 grating.

The method returns an astropy.table.Table
of the meta data.  If there were entries
from multiple :doc:`groups`, then the Table
is likely to be masked.  This table also
includes a new column 'Group' specifying
the group origin for each entry.

Here is another example::

    qdict = {'R': (4000.,1e9), 'WV_MIN': (0., 4000.)}
    qmeta2 = sdb.query_meta(qdict)

Now we are restricting on the spectral resolution
and wavelength coverage.

Meta Data from a Position on the Sky
------------------------------------

One can query the database for spectra
around a given location on the sky.
For convenience, the formatting of the
sky position includes many options.  One also inputs
a search radius which is either an Angle or Quantity.

Here is a simple example with a small search radius::

    meta = sdb.meta_from_position((0.0019,17.7737), 1*u.arcsec)

As with query_meta() from above, the meta_from_position()
method returns an astropy.table.Table with each row giving
the meta data of each spectrum matching the query.

One may commonly wish to restrict the query by data
:doc:`groups`.  Simply provide the list of groups::

    meta = sdb.meta_from_position((2.813500,14.767200), 20*u.deg, groups=['GGG','HD-LLS_DR1'])

The default is to return all spectra satisfying the position
query in each of the groups.  Set in_all_groups=True to require
that the source occur in all of the input groups.

Querying with a List of Coordinates
-----------------------------------

One can query the database with a set of coordinates,
each of which is matched to a small tolerance
(default: 0.5 arcseconds).
The input is an astropy.coordinate.SkyCoord object.
Here is an example::

    coords = SkyCoord(ra=[0.0028,0.0019], dec=[14.9747,17.7737], unit='deg')
    matches, subcat, IDs = sdb.qcat.query_coords(coords)

The outputs have the same size as the input set of coordinates
and are aligned.  As in the other queries, these are a bool array
indicating a match, the sub-catalog with rows ordered by the
input coordinates (non-matches are blank), and the IDKEY values.
Sources that do not match by coordinate have IDKEY=-1 and those
that match coordinates but fail some other criterion have
IDKEY=-2.

Here are a few other examples::

    qdict = dict(zem=(1.0, 2.5))
    matches, subcat, IDs = sdb.qcat.query_coords(coords, query_dict=qdict)

and::

    matches, subcat, IDs = sdb.qcat.query_coords(coords, groups=['BOSS_DR12'])


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
