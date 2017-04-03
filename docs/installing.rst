.. highlight:: rest

*****************
Installing specdb
*****************

This document describes how to install the `specdb`
repository.  We also describe
:ref:`download-public`.

Installing Dependencies
=======================
We have and will continue to keep the number of dependencies low.
There are a few standard packages that must be installed
and one package `linetools` under review for
`astropy` affiliated status.

In general, we recommend that you use Anaconda for the majority of
these installations.

Detailed installation instructions are presented below:

Python Dependencies
-------------------

specdb depends on the following list of Python packages.

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_
to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.12 or later
* `astropy <http://www.astropy.org/>`_ version 1.3 or later
* `scipy <http://www.scipy.org/>`_ version 0.18 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT5 <https://wiki.python.org/moin/PyQt/>`_ version 5 (needed for linetools)
* `h5py <https://www.h5py.org/>`_ version 2.6 (for data I/O)

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python|numpy|astropy|scipy|matplotlib|pyqt|h5py"

If the packages have been installed, this command should print
out all the packages and their version numbers.

If any of these packages are missing you can install them
with a command like::

	conda install h5py

If any of the packages are out of date, they can be updated
with a command like::

	conda update scipy

Installing linetools
--------------------
The latest version of `linetools <https://github.com/linetools/linetools/>`_
is also required for `specdb`. `linetools` is a package designed for the
analysis of 1-D spectra. The installation steps for `linetools` are
provided `here <http://linetools.readthedocs.io/en/latest/install.html/>`_.

Installing specdb
=================

Presently, you must download the code from github::

	#go to the directory where you would like to install specdb.
	git clone https://github.com/specdb/specdb.git

From there, you can build and install with

	cd specdb
	python setup.py install  # or use develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)


.. _download-public:

Downloading the public Data Bases
=================================

*specdb* provides a set of scripts to
:ref:`download-scipts`.  We summarize
the primary ones here.

igmspec
-------

`igmspec` is a public database intended to contain all published spectra associated
to intergalactic medium (IGM) studies (see
`paper_coming_soon <https://www.arxiv/>`_
for further details).
To grab `igmspec` you can download the database
(in HDF5 format) from its public
domain using the script provided by `specdb`.
You can download the file to
any location of your choice but we recommend to locate
it in the ./data/DB/ subdirectory within specdb, e.g.::

    cd specdb/data/DB
    specdb_get_igmspec

This will start a download of the most recent `igmspec` database; the current file is v02
and has ~26Gb. Once the file is downloaded,
you need to make sure the shell environmental
variable `$SPECDB` points to the host directory.

See :ref:`download-igmspec` for additional details on the script.
