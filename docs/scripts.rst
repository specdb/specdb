.. highlight:: rest

**************
specdb Scripts
**************

This file summarizes the specdb scripts
(used outside Python).  These are installed
within your standard Python script path (e.g.
~/anaconda/bin).

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Simple Scripts <Simple_Scripts>

plot_specdb
===========

Plot a spectrum at the given coordinate.  One can
restrict the database and/or surveys used and/or select the desired
spectrum from the available list.  By default, the
XSpecGui gui from linetools is called to display
the spectrum.   Here is the help::

   plot_specdb -h
    usage: plot_specdb [-h] [--tol TOL] [--meta] [-s SURVEY] [--select SELECT]
                       [--mplot MPLOT] [--db_file DB_FILE]
                       coord dbase

    plot_specdb script v0.3

    positional arguments:
      coord                 Coordinates, e.g. J081240.7+320809
      dbase                 Database [igmspec,all,priv]

    optional arguments:
      -h, --help            show this help message and exit
      --tol TOL             Maximum offset in arcsec [default=5.]
      --meta                Show meta data? [default: True]
      -s SURVEY, --survey SURVEY
                            Name of Survey to use
      --select SELECT       Index of spectrum to plot (when multiple exist)
      --mplot MPLOT         Use simple matplotlib plot [default: False]
      --db_file DB_FILE     Full path of db_file

Here is an example or two::

   plot_specdb J220248.31+123656.3 priv --db_file=qpq_optical.hdf5
   plot_specdb J220248.31+123656.3 igmspec


sdss_spec
=========

Grab data from the SDSS/BOSS survey with plate-fiber notation.
Here is the help::

   $sdss_spec -h
    usage: sdss_spec [-h] [-s SURVEY] [--select SELECT] [-p] plate fiberid dbase

    sdss_spec script v0.1

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

