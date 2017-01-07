.. highlight:: rest

*******
Spectra
*******

Description
===========

The spectra are stored in a series of HDF5 :doc:`groups`
as multi-dimension arrays (npixels x nspectra x nkeys).
Each array has the following keys.

=============  ======= =============================================
Key            Type    Description
=============  ======= =============================================
wave           float64 Wavelength array; default is Angstroms
flux           float32 Flux array; default is unitless
sig            float32 Error array; same units as flux
co (optional)  float32 Continuum array; same units as flux
=============  ======= =============================================

.. _XSpectrum1D: http://linetools.readthedocs.io/en/latest/xspectrum1d.html

The Python software included with `specdb` can
read these arrays into the
`XSpectrum1D`_ class provided
in `linetools <http://linetools.readthedocs.io/en/latest/>`_.

.. _retrieve-spectra:

Retrieving Spectra
==================

There are several methods to retrieve spectra
from a database instantiated in the :ref:`specdb-class`.
But we *strongly recommend* that the default approach be:

1. Perform a meta query.
2. Inspect the table returned.
3. Retrieve the spectra with this table (or a cut-down version).

A series of simple examples are given in the
`Retrieving Spectra Notebook <https://github.com/specdb/specdb/blob/master/docs/nb/Retrieve_Spectra.ipynb>`_.
The main methods are also described below, with examples.
The lower level methods using the :ref:`interface-group`
class are described in :ref:`group-spectra`.

Spectra from a Meta data Query
------------------------------

As noted above, this is the recommended approach to retrieving
spectra.  One is required to first generate a meta data Table
by :ref:`query-meta`.  One may then retrieve all of the
associated spectra.::

    # Meta query
    qdict = {'TELESCOPE': 'Gemini-North', 'NPIX': (1580,1583), 'DISPERSER': ['B600', 'R400']}
    qmeta = sdb.query_meta(qdict)

    # Retrieve spectra
    spectra = sdb.spectra_from_meta(qmeta)

This returns an `XSpectrum1D`_ object containing all of the
spectra associated with the meta data Table, aligned to the
Table.  Of course, one can first slice the meta data Table to
retrieve only a subset of the spectra.

Spectra for a single Source
---------------------------

A common usage of specdb may be to grab all of the spectra
related to a single source.  The method `spectra_from_coord`
takes an input coordinate (in a range of :ref:`coord_formats`),
identifies the closest catalog source within a given tolerance
(default is 0.5") and returns all of the spectra and meta data
within the database for that source.  Here is an example call::

   spec, meta = sdb.spectra_from_coord('J223438.52+005730.0')

spec and meta are a `XSpectrum1D`_ object and an
astropy.Table object.

One can restrict the call to grab spectra from a subset of the
groups or by querying the meta data, e.g.::

    qdict = dict(DISPERSER='R400')
    spec, meta = sdb.spectra_from_coord('001115.23+144601.8', query_dict=qdict)

In this example, we require the spectra have been taken with the R400
grating.  Note that queries are performed on both the source catalog
and the :doc:`meta` table of each group.

A related method is to retrieve all the spectra for a source
given its IDKEY::

    spec, meta = sdb.spectra_from_ID(3244)

The objects returned are the same as above.

Spectra from Coordinates in a Group
-----------------------------------

Another common usage may be to retrieve the spectra from a list of coordinates
from a single group.  The `spectra_in_group` method accomplishes this most
efficiently.  Here the input must be a SkyCoord object containing the
coordinates for one or more sources.  An example call::

    coords = SkyCoord(ra=[0.0028, 0.0019], dec=[14.9747, 17.77374], unit='deg')
    spec, meta= igmsp.spectra_in_group(coords, 'BOSS_DR12')

The output is an `XSpectrum1D`_ object containing the spectra and
an astropy.Table of the meta data.  This
returns only the first spectrum and meta row identified
in the group for each source, ordered the same as the input coordinates.

For cases where one or more spectra may be present, you may
wish to restrict by providing a :doc:`query_dict`, e.g.::

    coords = SkyCoord(ra=[2.8135,16.5802], dec=[14.7672, 0.8065], unit='deg')
    qdict = dict(DISPERSER='R400')
    spec, meta = sdb.spectra_in_group(coords, 'GGG', query_dict=qdict)

This requires the spectra returned were taken with the R400 grating.


**Note:** This method will raise an *IOError* if one or more of the input
coordinates are not within the requested group to within
the tolerance parameter (default = 0.5") or if one or more sources
fails to match an input :doc:`query_dict`.
