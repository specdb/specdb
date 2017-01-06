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

.. _group-query-meta:

Querying a Single Meta Table
============================

Using the :ref:`interface-group` class, one
queries a single meta data table using
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


