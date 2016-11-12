.. highlight:: rest

****************
Private Database
****************

It is possible within `specdb` to generate a private
database that can be used in tandem with the public
database(s). The following describes the procedure.

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Private <Private_Ingest>

Setup
=====

Directory Tree
--------------

The main input is a simple directory tree containing the
FITS files of individual spectra.  Each branch off the main
tree generates a unique dataset in the database.  It is also
expected (although not strictly required) that each branch
contains FITS files from a single instrument.  One is
allowed to have sub-folders in a branch, although this is
also not recommended.

Here is an example of a directory tree
(from the test dataset in specdb)::

   ├── privateDB
   |  ├── testDB_ztbl.fits
   |  ├── ESI
   |  |  ├── ESI_meta.json
   |  |  ├── SDSSJ220758.30+125944.3_F.fits
   |  |  ├── SDSSJ220758.30+125944.3_E.fits
   |  ├── LRIS
   |  |  ├── LRIS_meta.json
   |  |  ├── SDSSJ230044.36+015541.7_r600_F.fits
   |  |  ├── SDSSJ230044.36+015541.7_b400_F.fits
   |  ├── COS
   |  |  ├── COS_meta.json
   |  |  ├── J095240.17+515250.03.fits.gz
   |  |  ├── J095243.05+515121.15.fits.gz


Meta parameter file
-------------------

Every dataset in specdb must include a meta table and there
is a required set of columns (see :doc:`meta`) for the list
and descriptions.

The code can automatically generated the meta table from the
data files, but it highly recommended that you guide this process
by providing a meta parameter file in each branch of the tree.
It must be a JSON file and it must end in _meta.json.

Here is an example file::

   {
      "maxpix": 60000,
      "meta_dict": {
         "TELESCOPE": "HST"
      },
      "parse_head": {
         "DATE-OBS": "DATE",
         "GRATING": "OPT_ELEM",
         "INSTR": "INSTRUME",
         "R": true
      }
   }

This example sets the maximum number of pixels in any given file
within the branch (default: 10000).  The values in meta_dict are
set for each file.  The items in parse_head indicate which header
keyword to use for each meta parameter.

The spectral resolution may be dynamically calculated.  That
is the default and it is asserted with True in this example.

Redshift table
--------------

Each unique source in your database (RA, DEC) is required to also
have a redshift.   This must be supplied as a separate Table with
at least the following columns:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
RA          float    RA in degrees
DEC         float    DEC in degrees
ZEM         float    Redshift value
ZEM_SOURCE  str      Name of the source (e.g. SDSS, BOSS)
==========  ======== ============================================

One can generate one's own or specify any of the public specdb databases.
If you generate your own, place in the top-level of the tree and
give it an extension _ztbl.fits.

Quick go
========

The database construction is intended to be run in one go with
a single command from the command line.  One uses the specdb_privatedb
script.  Here is the current usage::


And here is an example of running it on the test DB::



Grab Files
----------

Meta
----

From the list of FITS files, a META table is generated.
This includes redshifts taken from the Myers catalog (when available)::

   meta = pbuild.mk_meta(spec_files, fname=True, skip_badz=True)

The *fname* flag indicates that the RA/DEC are to be parsed
from the FITS filename.  The *skip_badz* flag allows the code
to skip sources that are not cross-matched to the Myers catalog.

Spectra
=======

Spectra are simply ingested into an HDF5 file
provided they can be read with
`linetools.spectra.io.readspec`::

   pbuild.ingest_spectra(hdf, 'test', meta)

One Step
========

It is recommended that all of the above steps be run in
one go with the mk_db method::

   pbuild.mk_db([tree], ['test'], 'tmp.hdf5', ztbl, fname=True, skip_badz=True, nmax_pix=50000)

One inputs a list of directory trees and a list of names
for each one. Key words are passed to the various methods.
