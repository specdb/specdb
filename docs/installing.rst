.. highlight:: rest

*****************
Installing specdb
*****************

This document describes how to install specdb.

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

specdb depends on the following list of Python packages.

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_
to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.11 or later
* `astropy <http://www.astropy.org/>`_ version 1.1 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt/>`_ version 4 (needed for linetools)
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
is also required for specdb. linetools is a package designed for the
analysis of 1-D spectra. The installation steps for linetools are
provided `here <http://linetools.readthedocs.io/en/latest/install.html/>`_.

Installing specdb
=================

Presently, you must grab the code from github::

	#go to the directory where you would like to install PYPIT.
	git clone https://github.com/specdb/specdb.git

From there, you can build and install either with install or develop
(we recommend the latter), e.g.::

	cd specdb
	python setup.py develop


This should install the package and scripts.
Make sure that your PATH includes the standard
location for Python scripts (e.g. ~/anaconda/bin)


