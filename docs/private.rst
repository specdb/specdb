.. highlight:: rest

****************
Private Database
****************

It is possible within `specdb` to generate a private
database that can be used in tandem with the public
database(s). The following document describes the procedure.

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
tree generates a unique group in the database.  It is also
expected (although not required) that each branch
contains FITS files from a single instrument.  One is
allowed to have sub-folders in a branch, although this is
also not recommended.

Here is an example of a directory tree
(from the test dataset in specdb)::

   ├── test_privateDB
   |  ├── testDB_ztbl.fits
   |  ├── ESI
   |  |  ├── ESI_meta.json
   |  |  ├── ESI_meta.ascii
   |  |  ├── SDSSJ220758.30+125944.3_F.fits
   |  |  ├── SDSSJ220758.30+125944.3_E.fits
   |  ├── LRIS
   |  |  ├── LRIS_meta.json
   |  |  ├── SDSSJ230044.36+015541.7_r600_F.fits
   |  |  ├── SDSSJ230044.36+015541.7_b400_F.fits
   |  ├── COS
   |  |  ├── COS_meta.json
   |  |  ├── COS_ssa.json
   |  |  ├── J095240.17+515250.03.fits.gz
   |  |  ├── J095243.05+515121.15.fits.gz


Meta parameter file
-------------------

Every group in `specdb` includes a meta table and there
is a required set of columns (see :doc:`meta`) for the list
and descriptions.

The code can automatically generate the meta table from the
data files alone, but it is highly recommended that you
guide this process by providing a meta parameter file in each branch of the tree.
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

Meta table
----------

It is *recommended* that one attempt to extract as much of the meta
data from the spectra files as possible.  However, it is possible to
include additional meta data (or to over-ride the meta data in the
spectra) by including a meta data table.  The format is either
ASCII or FITS with a .ascii or .fits extension that must be
read by astropy.table.Table.read().

SPEC_FILE is a required column which gives the name of the spectral
file to match against the meta data.

Here is an example from the test suite (ESI_meta.ascii)::

   SPEC_FILE                         tGRB
   SDSSJ172524.66+303803.9_F.fits    2009-11-23:10:12:13.2
   SDSSJ220758.30+125944.3_F.fits    2007-08-13:10:22:23.3

This will add the time of the GRB to the meta data table.

SSA info
--------

*specdb* includes software to enable SSA queries of your
database.  For this to work, however, one must provide
a few additional fields for each data group.  These are
provided with a JSON file in each branch with extension
_ssa.json.

The required keys are:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
Title       str      Title for the data group
flux        str      Sets units and ucd for the flux.  Allowed values are
                     flambda, normalized
fxcalib     str      Sets Calibration field.  Allowed values are
                     NORMALIZED, ABSOLUTE, RELATIVE
==========  ======== ============================================

See the COS_ssa.json file in the test suite for an example.

One is also required to include a Publisher value.  This is defaulted
to 'Unknown', but can be set in the call to mk_db() or with the
--publisher keyword in the script.

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
                           [--version VERSION] [--fname]
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
      --fname            Parse RA/DEC from filename?

And here is an example of running it on the test DB::

   cd specdb/specdb/data/
   specdb_privatedb testDB test_privateDB tst_DB.hdf5

This will create a private DB called `testDB` from the directory tree
`test_privateDB`; the database itself is contained in a single .hdf5 named
`tst_DB.hdf5`

Within Python
-------------

Here is a call for the test database in one go
from within Python::

   from specdb.build import privatedb as pbuild
   # Read z table
   ztbl = Table.read(specdb.__path__[0]+'/data/test_privateDB/testDB_ztbl.fits')
   # Go
   tree2 = specdb.__path__[0]+'/data/test_privateDB'
   pbuild.mk_db(tree2, 'testDB', 'tst_DB.hdf5', ztbl, fname=True)

If ztbl == 'igmspec',  the code will attempt to load the
IgmSpec database and use the quasars catalog for redshifts.

Step by Step in Python
======================

It is possible that you will need to customize things
further.  This section describes the step-by-step
approach from within Python.
Note that there is additional code within mk_db() that may
be required.

Get Started
-----------

Start the main catalog and set your private ID_KEY::

   id_key = 'TEST_ID'
   maindb, tkeys = spbu.start_maindb(id_key)

This sets a global variable within specdb.build.utils

Grab Files
----------

The grab_files() method searches through a given branch to find
all FITS files and a meta parameter file.  By default, the code
ignores any files with the following extensions:
'c.fits', 'C.fits', 'e.fits', 'E.fits', 'N.fits', 'old.fits'.

Here is an example call::

   branch = specdb.__path__[0]+'/data/test_privateDB/ESI/'
   flux_files, meta_file, custom_meta_table = pbuild.grab_files(branch)


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


Add Groups and IDs
------------------

Update the main catalog and fiddle about with a few tags in the
meta table.  Also update the group dict::

   gdict = {}
   flag_g = spbu.add_to_group_dict('COS', gdict)
   maindb = pbuild.add_ids(maindb, meta, flag_g, tkeys, id_key, first=(flag_g==1))

The group dict is eventually written to the HDF5 file.

Ingest Spectra
--------------

The ingest_spectra() method loops through the spectral files,
reads each, and populates a hdf5 dataset.  It also converts
the meta table into a separate, parallel hdf5 dataset.

Here is an example::

   hdf = h5py.File('tmp.hdf5','w')
   pbuild.ingest_spectra(hdf, 'test', meta, max_npix=50000)


Finish
------

Before writing, the code tests whether the meta data tables
can be stacked using the specdb.utils.clean_vstack() method.
This may be required for queries of the meta data and spectra
extraction.  The code will hit a pdb.set_trace() if this fails.

Write the catalog and close the HDF5 file.::

   zpri = [str('SDSS'), str('BOSS')]
   pbuild.write_hdf(hdf, 'TEST_DB', maindb, zpri, gdict, 'v01')

Also writes the creation date and sets the version.
*zpri* is a list of strings indicating the priority given to
redshift assignment.
