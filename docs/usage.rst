.. highlight:: rest

*****************
Usage with Python
*****************

This file introduces the :ref:`specdb-class` and provides
links to the main actions
for using a `specdb` database from within Python.

See the :doc:`scripts` documentation for a discussion of
command-line usage from the command line.

Overview
========

There are three standard activities related to using
a specdb database:

1. :ref:`query-catalog`

2. :ref:`query-meta`

3. :ref:`retrieve-spectra`

Use the above links to examine
the documentation and examples
for each of these activities.

.. _specdb-class:

SpecDB Class
============

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


