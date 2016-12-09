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

   ├── test_privateDB
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

Every dataset in `specdb` must include a meta table and there
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
also set for all the files in the branch. The items in
parse_head indicate which header keyword to use for retrieving
the corresponding values for each individual file of the branch. These
values may be different from file to file and this is a convenient way to get
them directly from the headers of the FITS files. The spectral resolution (`R`) may
be dynamically calculated within `specdb`.  That is the default and it is asserted
with `True` in this example.

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

One can generate one's own table or specify any of the public `specdb`
databases (e.g. `igmspec`). If you generate your own, place in the top-level
of the database tree and give it an extension _ztbl.fits (see example of
tree structure above).

Spectra
-------

Spectra will be ingested provided they can be read with
`linetools.spectra.io.readspec`.  You can test whether this is the case
by running::

   lt_xspec name_of_spectrum

on any of your files.  The default is to take any FITS file
in the branch (and sub-folders) except those files with these
extensions: 'c.fits', 'C.fits', 'e.fits', 'E.fits', 'N.fits', 'old.fits'.

Quick go
========

Script
------

The database construction is intended to be run in one go with
a single command from the command line. One uses the specdb_privatedb
script. Here is the current usage::

   usage: specdb_privatedb [-h] [--ztbl ZTBL] [--zspecdb ZSPECDB]
                        db_name tree_path outfile

   Generate a private specdb DB

   positional arguments:
      db_name            Name of your private DB
      tree_path          Path to the directory tree of spectral files
      outfile            Filename for the private DB HDF5

   optional arguments:
      -h, --help         show this help message and exit
      --ztbl ZTBL        Name of data file containing redshift info
      --zspecdb ZSPECDB  Name of specdb DB to use for redshifts
      --version VERSION  Version of the DB; default is `v00`

And here is an example of running it on the test DB::

   cd specdb/specdb/data/
   specdb_privatedb testDB test_privateDB tst_DB.hdf5

This will create a private DB called `testDB` from the directory tree
`test_privateDB`; the database itself is contained in a single .hdf5 named
`tst_DB.hdf5`

Within Python
-------------

It is possible you will need to customize things to the
extent that you will want to generate the database
from inside Python.

Here is a call for the test database::

   from specdb.build import privatedb as pbuild
   # Read z table
   ztbl = Table.read(specdb.__path__[0]+'/data/test_privateDB/testDB_ztbl.fits')
   # Go
   tree2 = specdb.__path__[0]+'/data/test_privateDB'
   pbuild.mk_db(tree2, 'testDB', 'tst_DB.hdf5', ztbl, fname=True)

If ztbl == 'igmspec',  the code will attempt to load the
IgmSpec database and use the quasars catalog for redshifts.

Main Steps
==========

Here are the main, individual steps taken when generating the
private DB.  There is, however, additional code within mk_db()
that is required.


Grab Files
----------

The grab_files() method searches through a given branch to find
all FITS files and a meta parameter file.  By default, the code
ignores any files with the following extensions:
'c.fits', 'C.fits', 'e.fits', 'E.fits', 'N.fits', 'old.fits'.

Here is an example call::

   branch = specdb.__path__[0]+'/data/test_privateDB/ESI/'
   flux_files, meta_file = pbuild.grab_files(branch)


Meta
----

From the list of FITS files, a META table is generated.
The redshift table must be supplied (as an astropy Table).
Here is an example call::

   meta = pbuild.mk_meta(flux_files, ztbl, fname=True, mdict=mdict, parse_head=pdict, skip_badz=True)

The *fname* flag indicates that the RA/DEC are to be parsed
from the FITS filename.  The *skip_badz* flag allows the code
to skip sources that are not cross-matched to redshift table
(instead of terminating).


Ingest Spectra
--------------

The ingest_spectra() method loops through the spectral files,
reads each, and populates a hdf5 dataset.  It also converts
the meta table into a separate, parallel hdf5 dataset.

Here is an example::

   hdf = h5py.File('tmp.hdf5','w')
   pbuild.ingest_spectra(hdf, 'test', meta, max_npix=50000)

