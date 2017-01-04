.. highlight:: rest

*********************
Basic Usage in Python
*********************

This file summarizes some of the simple cases
of using a `specdb` database from within Python.
See the :doc:`scripts` documentation for a discussion of
command-line usage from the command line.

Overview
========

There are three standard activities related to using
a specdb database:

1. :ref:`query-catalog`

2. Querying the meta tables

3. Retrieving spectra

We provide separate documentation and examples
for each of these activities.

SpecDB
======

Essentially all of database interfacing is performed
through the SpecDB Class (or one of its children).
This class loads the database, performs queries on
the catalogs and meta data, and retrieves the spectra.

Instantiation
-------------

Step 0 is to instantiate by pointing the Class
at the database of interest.  Here is an example::

    from specdb.specdb import SpecDB
    specdb = SpecDB(db_file='/my_path/IGMspec_DB_v02.hdf5')

At instantiation, the SpecDB object performs the actions:

 - Opens a pointer to the HDF5 file
 - Loads a QueryCatalog object into specdb.qcat
 - Loads the source catalog as an astropy.Table in specdb.cat
 - Loads the list of data groups in specdb.groups
 - Loads a *dict* translating the data group flags in specdb.group_dict
 - Generates an empty *dict* in specdb._gdict which will be used to interface with the data groups


Children
--------

The various public databases may have a unique child
of SpecDB.  Currently, there is:

========== ====================================================
Database   SpecDB Child
========== ====================================================
igmspec    IgmSpec
========== ====================================================

Here is an example instantiation for *igmspec*::

    igmsp = IgmSpec()  # Loads the highest version in $IGMSPEC_DB

This loads the highest version number of any HDF5 files located
in the $IGMSPEC_DB folder.


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

