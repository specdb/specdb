.. highlight:: rest

*********************
Coordinates in specdb
*********************

The algorithms in specdb primarily make use of the SkyCoord
class in astropy.  For input in scripts and methods, however,
most of the routines take a range of formats which are converted
by linetools.utils.radec_to_coord to a SkyCoord instance.


.. _coord_formats:

Formats
=======

Here are the generally allowed formats for input
coordinates:

======= ========= ======================== =======================
Style   Type      Example                  Comment
======= ========= ======================== =======================
J2000   str       'J123411.12+123411.1'    The 'J' is not required
colon   str       '00:22:33.1,+12:22:33.3'
radec   tuple     (23.2311, -12.2311)      RA, DEC in deg
coord   SkyCoord
======= ========= ======================== =======================
