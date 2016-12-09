.. highlight:: rest

************
Using specdb
************

This file summarizes some of the simple usage cases
of `specdb` from within Python.
See the :doc:`scripts` documentation for a discussion of
command-line usage cases outside of Python.

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Simple Usage <Simple_Usage>

Setup
=====

The main class of specdb is simply instantiated::

    from specdb.specdb import SpecDB
    igmsp = SpecDB(db_file='/raid/IGMSPEC_DB/IGMspec_DB_v02.hdf5')

This loads the database catalog::

    igmsp.qcat
    <QueryCatalog:  DB_file=/u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 with 377018 sources Loaded groups are [u'BOSS_DR12', u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7'] >

and opens the hdf5 file (without loading the data).


Querying the Catalog
====================

There are several methods that interface with the primary
source catalog.

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
an ndarray of IDs with shape matching the input list.
Coordinates without a match within the tolerance
have -1 values .  Here is an example::

    coords = SkyCoord(ra=[0.0019,1.2321], dec=[17.7737,-12.2332], unit='deg')
    IDs = igmsp.qcat.match_coord(coords)

One can further restrict the search to a specific group

    IDs = igmsp.match_coord(coords, group='BOSS_DR12')

Sources that are a match in position but not within the group
have an ID=-2.

show_cat
--------

A printout of the catalog values for a list of IDs is provided
by `show_cat`::

   igmsp.qcat.show_cat(IDs)

This includes the flag_group values which indicate the groups
that include a given source.  The catalog only shows a single
entry per source and only those sources with ID values within
the catalog (e.g. negative values are ignored).

Grabbing Spectra
================

allspec_at_coord
----------------

A common usage of specdb may be to grab all of the spectra
related to a single source.  The method `allspec_at_coord`
takes an input coordinate (in a range of :ref:`coord_formats`),
identifies the closest catalog source within a given tolerance
(default is 0.5") and returns all of the spectra and meta data
within the database for that source.  Here is an example call::

   speclist, metalist = igmsp.allspec_at_coord('J223438.52+005730.0')

speclist and metalist are lists of XSpectrum1D and astropy.Table objects,
one for each group that includes the source.

One can restrict the call to grab spectra from a subset of the
groups, e.g.::

   speclist, metalist = igmsp.allspec_at_coord('J223438.52+005730.0', igroup=['HD-LLS_DR1'])
   spec = speclist[0]

coords_to_spec
--------------

Another common usage will be to grab the spectra for a list of coordinates
from a single group.  The `coords_to_spec` method accomplishes this most
efficiently.  Here the input must be a SkyCoord object containing the
coordiantes for one or more sources.  An example call::

    coords = SkyCoord(ra=[0.0028, 0.0019], dec=[14.9747, 17.77374], unit='deg')
    spec, meta= igmsp.coords_to_spectra(coords, 'BOSS_DR12')

The output is an XSpectrum1D object containing the spectra and
an astropy.Table of the meta data.  The default mode is to
return the first spectrum and meta row in the group for each
source, ordered the same as the input coordinates.

Alternatively, you may request all of the spectra matching to the
input set of IDs with `all_spec=True`.  In this case, the ordering
is simply how the group data were ingested.

**Note:** This method will raise an error if one or more of the input
coordinates are not within the requested group to within
the tolerance parameter (default = 0.5").

