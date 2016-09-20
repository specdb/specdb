.. highlight:: rest

*************
Using igmspec
*************

This file summarizes some of the simple usage cases
of igmspec from within Python.
See the :doc:`scripts` documentation for a discussion of
command-line usage cases outside of Python.

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Simple Usage <Simple_Usage>

Setup
=====

The main class of igmspec is simply instantiated::

    from igmspec.igmspec import IgmSpec
    igmsp = IgmSpec()

This loads the database catalog::

    igmsp.qcat
    <QueryCatalog:  DB_file=/u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 with 377018 sources Loaded surveys are [u'BOSS_DR12', u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7'] >

and opens the hdf5 file (without loading the data)::

    igmsp.idb
    <InterfaceDB:  DB_file=/u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 Loaded surveys are [u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7']>


Grabbing Spectra
================

A common usage of igmspec may be to grab a single spectrum from
a single survey.  Here is an example::

   speclist, metalist = igmsp.spec_from_coord('J223438.52+005730.0', isurvey=['HD-LLS_DR1'])
   spec = speclist[0]


