.. highlight:: catalog

*******
Catalog
*******

Fundamental to a specdb database is the source catalog, written
to the 'catalog' dataset in the hdf5 file.  At a minimum this
table contains the RA, DEC, redshift info, a group flag,
and a catalog IDKEY.


Table
=====

The table comprising the source catalog has the following entries:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IDKEY       int      Unique identifier;  name of the key should include ID (e.g. IGM_ID)
flag_group  int      Bitwise flag indicating the groups that the source has spectra in
zem         float    Emission redshift of background source
sig_zem     float    Estimated error in the redshift
flag_zem    str      Key indicating source of the redshift (e.g. BOSS_PCA)
RA          float    Right Ascension (deg)
DEC         float    Declination (deg)
STYPE       str      Type of the source (e.g. QSO)
==========  ======== ============================================


