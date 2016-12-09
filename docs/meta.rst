.. highlight:: rest

*********
Meta Data
*********

Each group in specdb includes a meta data table
with each row describing each ingested spectrum.

Here are the required columns:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IDKEY       int      Database ID values.
zem         float    Emission redshift of background source
RA          float    Right Ascension (deg)
DEC         float    Declination (deg)
EPOCH       float    Coordinate epoch (typically 2000)
DATE-OBS    str      Date observed (YYYY-MM-DD)
R           float    Instrument resolution, :math:`\lambda/\Delta\lambda` (FWHM)
WV_MIN      float    Minimum wavelength of the spectrum
WV_MAX      float    Maximum wavelength of the spectrum
NPIX        int      Number of pixels in the spectrum; may include null values
GROUP_ID    int      Unique identifier for each row in the meta table
SPEC_FILE   str      Spectrum file name
INSTR       str      Instrument file name (see `Instruments and Gratings`_ below for definitions)
GRATING     str      Grating name (see `Instruments and Gratings`_ below for definitions)
TELESCOPE   str      Telescope name (see `Telescopes`_ below for definitions)
==========  ======== ============================================

Each specdb database should use a unique ID key for the Database ID
values.  For example, igmspec uses IGM_ID.

Instruments and Gratings
------------------------

The instruments used in specdb are provided in specdb.defs.instruments.
The following Table summarizes and defines the instruments
currently used within the specdb software:

==========  ========     ============================================
Instrument  Gratings     Description
==========  ========     ============================================
ACS         PR200L       HST/ACS Prism mode (slitless)
BOSS        BLUE         Blue channel spectrograph
 ..         RED          Red channel spectrograph
 ..         BOTH         Spectrum includes data from both spectrographs
COS         G130M        HST/COS spectrometer
 ..         G160M        HST/COS spectrometer
 ..         G130M/G160M  HST/COS spectrometer with combined gratings
 ..         G130M-G160M  HST/COS spectrometer with combined gratings
ESI         ECH          Echelette mode on Keck/ESI instrument
FOS         G160L        HST/FOS spectrometer
 ..         G130H        ..
 ..         G190H        ..
 ..         G270H        ..
FUSE        LWRS         FUSE spectrometer
GMOS-N      R400         Gemini North GMOS spectrometer
 ..         B600         ..
GMOS-S      R400         Gemini South GMOS spectrometer
 ..         B600         ..
GNIRS       ECH          Gemini GNIRS spectrometer
HIRES       BLUE         Blue cross-disperser on HIRES (aka HIRESb)
 ..         UV           Blue cross-dispereser on HIRES (historic name)
 ..         RED          Red cross-dispereser on HIRES (aka HIRESr)
 ..         BOTH         Spectrum includes data from both cross-dispersers
ISAAC       SW_MRes      VLT/Isaac spectrometer
LRISb       400/3400     Keck/LRIS blue camera
 ..         600/4000     ..
 ..         1200/3400    ..
LRISr       600/7500     Keck/LRIS red camera
 ..         400/8500     ..
 ..         1200/7500    ..
MagE        N/A          MagE spectrometer
MIKEb       BLUE         Blue camera of MIKE spectrometer
MIKE-Blue   BLUE         Alternative instrument name
MIKEr       RED          Red camera of MIKE spectrometer
MIKE-Red    RED          Alternative instrument name for MIKE
MIKE        BOTH         Spectrum is a splice of MIKEb and MIKEr data
MMT         ??           Defaults to MMT red channel spectrograph (RCS)
mmtbluechan 500GPM       MMT blue channel spectrograph (BCS)
NIRI        Hgrism_G5203 Gemini NIRI spectrometer
 ..         Kgrism_G5204 ..
NIRSPEC     Low-Res      Low resolution mode of NIRSPEC
MODS1B      G400L        MODS spectrometer on LBT; blue side
MODS1R      G670L        MODS spectrometer on LBT; red side
MOSFIRE     H            H-band mode
SDSS        BLUE         Blue channel spectrograph
 ..         RED          Red channel spectrograph
 ..         BOTH         Spectrum includes data from both spectrographs
STIS        G140L        HST/STIS spectrometer
 ..         G230L        ..
 ..          ..          And many more..
TSPEC       ECH          Palomar Triplespec spectrometer
XSHOOTER    UVB          VLT/XShooter spectrometer in UVB camera
 ..         VIS          VLT/XShooter spectrometer in VIS camera
 ..         NIR          VLT/XShooter spectrometer in NIR camera
UVES        BOTH         VLT/UVES spectrometer
WFC3        G280         WFC3 grism (slitless)
2dF         300B         Blue channel spectrograph
==========  ========     ============================================

Telescopes
----------

Here are the telescopes currently incorporated in specdb:

==============  ====================================================
Telescope       Website
==============  ====================================================
Gemini-N        http://www.gemini.edu
Gemini-S        http://www.gemini.edu
HST             http://www.stsci.edu/hst/
Keck-I          http://www.keckobservatory.org/
Keck-II         http://www.keckobservatory.org/
LBT             http://www.lbto.org/
Magellan/Clay   http://obs.carnegiescience.edu/Magellan
Magellan/Baade  http://obs.carnegiescience.edu/Magellan
MMT             https://www.mmto.org/
SDSS 2.5-M      https://www.sdss3.org/instruments/telescope.php
UKST            https://www.aao.gov.au/about-us/uk-schmidt-telescope
VLT             http://www.eso.org/public/teles-instr/paranal/
==============  ====================================================

