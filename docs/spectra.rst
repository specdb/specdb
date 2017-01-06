.. highlight:: rest

*******
Spectra
*******

The spectra are stored in a series of hdf5 datasets
as multi-dimension arrays.  This array has the following
format.

=============  ======= =============================================
Key            Type    Description
=============  ======= =============================================
wave           float64 Wavelength array; default is Angstroms
flux           float32 Flux array; default is unitless
sig            float32 Error array; same units as flux
co (optional)  float32 Continuum array; same units as flux
=============  ======= =============================================

The software included with specdb read these arrays into
a XSpectrum1D object from
`linetools <http://linetools.readthedocs.io/en/latest/>`_.


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
