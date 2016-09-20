.. highlight:: rest

****************
Private Database
****************

It is possible within igmspec to generate a private
database that can be used in tandem with the public
database.  The following describes the procedure.

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Private <Private_Ingest>

Setup
=====

The main input is a directory tree containing the
FITS files of individual spectra.  It may contain
sub-folders, e.g.::

    tree = os.getenv('DROPBOX_DIR')+'/QSOPairs/data/MMT_redux/'

From this tree, a list of flux files is generated::

    from igmspec.privatedb import build as pbuild
    flux_files = pbuild.grab_files(tree)

Meta
====

From the list of FITS files, a META table is generated.
This includes redshifts taken from the Myers catalog.::

   meta = pbuild.mk_meta(flux_files, fname=True, skip_badz=True)

The *fname* flag indicates that the RA/DEC are to be parsed
from the FITS file.  The *skip_badz* flag allows the code
to skip sources that are not cross-matched to the Myers catalog.

Spectra
=======

Spectra are simply ingested into an HDF5 file
provided they can be read with
linetools.spectra.io.readspec::

   pbuild.ingest_spectra(hdf, 'test', meta)

One Step
========

It is recommended that all of the above steps be run in
one go with the mk_db method::

   pbuild.mk_db([tree], ['test'], 'tmp.hdf5',fname=True, skip_badz=True)

One inputs a list of directory trees and a list of names
for each one.  Key words are passed to the various methods.

