.. highlight:: groups

******
Groups
******

This document describes the meta data of
a `specdb` database and methods to query it.

.. _meta-desc:

Description
===========

Each group in a `specdb` database contains a meta data
table and the corresponding set of spectra.
The two are aligned, row by row.

Within the HDF5 file, each group is stored
in a data group which contains two datasets:

1. meta data in hdf[group+'/meta']
2. spectra in hdf[group+'/spec']

.. _interface-group:

InterfaceGroup
==============

The InterfaceGroup Class is used to interface
to a group.  It stores the meta data, provides
a method to query the data, and provides methods
to extract the spectra.

The class is instantiated with the point to an
HDF5 object, the name of the datagroup, and
the IDKEY::

    igroup = InterfaceGroup(hdf, 'GGG', 'IGM_ID')

Upon instantiation, the meta data table is stored
in the .meta attribute as an astropy.table.Table.

.. _query-meta:

Querying a Meta Table
=====================

Using the :ref:`interface-group` class, on
queries a meta data table with
a :doc:`query_dict` object.
Here is an example::

    qdict = {'TELESCOPE': 'Gemini-North', 'INSTR': 'GMOS-N', 'NPIX': (1580,1583), 'DISPERSER': ['B600', 'R400']}
    matches, sub_meta, IDs = igroup.query_meta(qdict)

This query requests all spectra taken at
the Gemini-North telescope with the GMOS-N instrument
having between 1580-1583 pixels (inclusive)
and with either the B600 or R400 grating.

The method returns a bool array (matches) indicating which
rows of the catalog are matched, the sub-table of those rows,
and the IDKEY values of those rows.

Querying a Position on the Sky
------------------------------

One can query the database around a given location
on the sky.  For convenience, the formatting of the
sky position includes many options.  One also inputs
a search radius which is either an Angle or Quantity
from astropy.

Here is a simple example with a small search radius::

    matches, sub_cat, IDs = sdb.qcat.query_position('001115.24+144601.9', 10*u.arcsec)

The objects returned are a bool array (matches) indicating
which rows of the catalog matched, the sub-catalog of
those rows ordered by separation from the search position,
and the IDs of those sources also ordered by search position.

Here is an example with a wider search and restricting to
sources that have spectra in at least one of a set of groups::

    matches, sub_cat, IDs = sdb.qcat.query_position((2.813500,14.767200), 20*u.deg, groups=['SDSS_DR7','GGG'])

Here the input was an (ra,dec) tuple assumed to be in decimal degrees.
Finally, an example that includes a :doc:`query_dict` to further
refine the search (on emission redshift)::

    qdict = dict(zem=(1.0, 3.))
    matches, sub_cat, IDs = sdb.qcat.query_position('001115.24+144601.9', 20*u.deg, query_dict=qdict)

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
