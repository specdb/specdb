.. highlight:: rest

**************
specdb Scripts
**************

This file summarizes the *specdb* scripts
which are executed outside a Python shell,
but use Python code.
These are installed
within your standard Python script path (e.g.
~/anaconda/bin).

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Simple Scripts <Simple_Scripts>

.. _download-scripts:

Download a Database
===================

*specdb* includes a set of scripts for downloading its
public databases.  Here is a brief description of each.

.. _download-igmspec:

specdb_get_igmspec
------------------

The download script for *igmspec* is named specdb_get_igmspec.
Here is its usage::

    wolverine-6.local> specdb_get_igmspec -h
    usage: specdb_get_igmspec [-h] [-v VERSION]

    Grab the IGMspec DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate

The script uses a simple wget command to the URL of the
*igmspec* database.  It defaults to the most recent
version of the database.

specdb_get_uvqs
---------------

The download script for *UVQS* is named specdb_get_uvqs.
Here is its usage::

    wolverine> specdb_get_uvqs -h
    usage: specdb_get_uvqs [-h] [-v VERSION]

    Grab the UVQS DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate

The script uses a simple wget command to the URL of the
*UVQS* database.  It defaults to the most recent
version of the database.

specdb_plot
===========

Plot a spectrum at the given coordinate.  One can
restrict the database and/or surveys used and/or select the desired
spectrum from the available list.  By default, the
XSpecGui gui from linetools is called to display
the spectrum.   Here is the usage::

   wolverine-6.local> specdb_plot -h
    usage: specdb_plot [-h] [--tol TOL] [--meta] [-g GROUP] [--select SELECT]
                       [--mplot MPLOT] [--db_file DB_FILE]
                       coord dbase

    specdb_plot script v0.3

    positional arguments:
      coord                 Coordinates, e.g. J081240.7+320809 or 122.223,-23.2322
                            or 07:45:00.47,34:17:31.1
      dbase                 Database [igmspec,all,priv]

    optional arguments:
      -h, --help            show this help message and exit
      --tol TOL             Maximum offset in arcsec [default=5.]
      --meta                Show meta data? [default: True]
      -g GROUP, --group GROUP
                            Name of Group to use (e.g. BOSS_DR12)
      --select SELECT       Index of spectrum to plot (when multiple exist)
      --mplot MPLOT         Use simple matplotlib plot [default: False]
      --db_file DB_FILE     Full path of db_file

Here is an example or two::

   specdb_plot J220248.31+123656.3 priv --db_file=qpq_optical.hdf5
   specdb_plot J220248.31+123656.3 igmspec


sdss_spec
=========

Grab data from the SDSS/BOSS survey with plate-fiber notation.
Here is the help::

   $specdb_sdss -h
    usage: specdb_sdss [-h] [-s SURVEY] [--select SELECT] [-p] plate fiberid dbase

    specdb_sdss script v0.1

    positional arguments:
      plate                 Plate
      fiberid               FiberID
      dbase                 Database [igmspec,all]

    optional arguments:
      -h, --help            show this help message and exit
      -s SURVEY, --survey SURVEY
                            Name of Survey to use (BOSS_DR12 or SDSS_DR7)
      --select SELECT       Index of spectrum to plot (when multiple exist)
      -p, --plot            Plot with lt_xspec

Here is an example::

   specdb_sdss 377 321 igmspec


