
Simple Scripts with SpecDB (v1.2.1)
===================================

:download:`Download <nb/Simple_Scripts.ipynb>` this notebook.

.. code:: python

    # imports

Downloading a DataBase
----------------------

After installing specdb, you can grab the latest (or any previous)
version of a given DB with its *get\_xxxx* script. Here is the call for
IgmSpec:

::

    usage: specdb_get_igmspec [-h] [-v VERSION]

    Grab the IGMspec DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate

IGMSpec :: specdb\_get\_igmspec
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

v02
^^^

::

    specdb_get_igmspec -v v02

--------------

Checking a DB file :: specdb\_chk
---------------------------------

Script to check the version and date of a given DB file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Usage
~~~~~

::

    wolverine-4.local> specdb_chk -h
    usage: specdb_chk [-h] db_file

    Check a specdb DB file

    positional arguments:
      db_file     Database file

    optional arguments:
      -h, --help  show this help message and exit

Example
~~~~~~~

::

    profx.ucolick.org> specdb_chk ~/local/Python/igmspec/DB/IGMspec_DB_v02.hdf5
    specdb DB file is from the igmspec database
    specdb DB file version=v02 was created on 2016-Oct-25
    Latest version for specdb DB type=igmspec is version=v02
    Latest creation date for this DB version was 2016-10-25
    Oldest valid DB file for this DB version was 2016-10-10
    Dataset: 2QZ
    Dataset: BOSS_DR12
    Dataset: COS-Dwarfs
    Dataset: COS-Halos
    Dataset: ESI_DLA
    Dataset: GGG
    Dataset: HD-LLS_DR1
    Dataset: HDLA100
    Dataset: HSTQSO
    Dataset: HST_z2
    Dataset: KODIAQ_DR1
    Dataset: MUSoDLA
    Dataset: SDSS_DR7
    Dataset: XQ-100
    Dataset: catalog
    Dataset: quasars

--------------

Plot :: specdb\_plot
--------------------

::

    usage: specdb_plot [-h] [--tol TOL] [--meta] [-s SURVEY] [--select SELECT]
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
                            Name of Group to use
      --select SELECT       Index of spectrum to plot (when multiple exist)
      --mplot MPLOT         Use simple matplotlib plot [default: False]
      --db_file DB_FILE     Full path of db_file

Examples
~~~~~~~~

FJ0812+32
^^^^^^^^^

::

    specdb_plot J081240.7+320808 igmspec --group KODIAQ_DR1

J001115.23+144601.8
^^^^^^^^^^^^^^^^^^^

::

    specdb_plot J001115.23+144601.8 igmspec

--------------

Interface with SDSS/BOSS Database :: specdb\_sdss
-------------------------------------------------

::

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

Example
~~~~~~~

::

    UNIX> specdb_sdss 434 555 igmspec -p

