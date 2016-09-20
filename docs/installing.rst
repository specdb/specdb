.. highlight:: rest

******************
Installing igmspec
******************

This document describes how to install igmspec and
its database (DB).

Installing Dependencies
=======================
We have and will continue to keep the number of dependencies low.
There are a few standard packages that must be installed
and one astropy affiliated (soon) package -- linetools.

In general, we recommend that you use Anaconda for the majority of
these installations.

Detailed installation instructions are presented below:

Python Dependencies
-------------------

igmspec depends on the following list of Python packages.

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_
to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.11 or later
* `astropy <http://www.astropy.org/>`_ version 1.1 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt/>`_ version 4 (needed for linetools)

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python|numpy|astropy|scipy|matplotlib|PyQT"

If the packages have been installed, this command should print
out all the packages and their version numbers.

If any of these packages are missing you can install them
with a command like::

	conda install PyQT

If any of the packages are out of date, they can be updated
with a command like::

	conda update scipy

Installing linetools
--------------------
The latest version of `linetools <https://github.com/linetools/linetools/>`_
is also required for igmspec. linetools is a package designed for the
analysis of 1-D spectra. The installation steps for linetools are
provided `here <http://linetools.readthedocs.io/en/latest/install.html/>`_.

Installing igmspec
==================

Presently, you must grab the code from github::

	#go to the directory where you would like to install PYPIT.
	git clone https://github.com/pyigm/igmspec.git

From there, you can build and install either with install or develop
(we recommend the latter), e.g.::

	cd igmspec
	python setup.py develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)

Downloading the DataBase (DB)
=============================

The igmspec package comes with a script for downloading the
latest version of the DB (or previous versions).  It is named
`get_igmspec`.  To grab the latest version, go to the directory
where you want it (especially if you already have a version) and do::

    get_igmspec

It will launch a wget command.  The size of iIGMspec_DB_v01.hdf5
is 5.8Gb.  You will then need to set the
IGMSPEC_DB environmental variable (if you have not already done so),
e.g.::

    setenv IGMSPEC_DB /u/xavier/local/Python/igmspec/DB/

Enjoy!
