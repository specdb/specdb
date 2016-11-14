.. highlight:: rest

************
Spectra Data
************

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
