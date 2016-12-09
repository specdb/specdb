
Ingesting Private Datasets (v1.1)
=================================

:download:`Download <nb/Private_Ingest.ipynb>` this notebook.

:download:`Download <nb/Private_Ingest.ipynb>` this notebook.

:download:`Download <examples/Private_Ingest.ipynb>` this notebook.

.. code:: python

    # imports
    import h5py
    from igmspec.privatedb import build as pbuild

Test on Single Folder
---------------------

.. code:: python

    tree = os.getenv('DROPBOX_DIR')+'/QSOPairs/data/MMT_redux/'

.. code:: python

    reload(pbuild)
    flux_files = pbuild.grab_files(tree)
    len(flux_files)




.. parsed-literal::

    87



.. code:: python

    flux_files[:5]




.. parsed-literal::

    [u'/u/xavier/Dropbox//QSOPairs/data/MMT_redux//./SDSSJ001605.89+005654.3_b800_F.fits',
     u'/u/xavier/Dropbox//QSOPairs/data/MMT_redux//./SDSSJ001607.27+005653.1_b800_F.fits',
     u'/u/xavier/Dropbox//QSOPairs/data/MMT_redux//./SDSSJ001743.07+050327.3_b500_F.fits',
     u'/u/xavier/Dropbox//QSOPairs/data/MMT_redux//./SDSSJ001744.13+050334.8_b500_F.fits',
     u'/u/xavier/Dropbox//QSOPairs/data/MMT_redux//./SDSSJ001807.36+161257.5_b500_F.fits']



Directory Tree
--------------

.. code:: python

    tree2 = os.getenv('DROPBOX_DIR')+'/QSOPairs/data/'

.. code:: python

    reload(pbuild)
    mflux_files = pbuild.grab_files(tree2)
    len(mflux_files)




.. parsed-literal::

    25451



.. code:: python

    mflux_files[:5]




.. parsed-literal::

    [u'/u/xavier/Dropbox//QSOPairs/data/BOSS/BOSSLyaDR9_cat.fits.gz',
     u'/u/xavier/Dropbox//QSOPairs/data/BOSS/BOSSLyaDR9_cat.fits.gz',
     u'/u/xavier/Dropbox//QSOPairs/data/ESI_redux/SDSSJ001813.89+142455.6_F.fits.gz',
     u'/u/xavier/Dropbox//QSOPairs/data/ESI_redux/SDSSJ002801.18-104933.9_F.fits.gz',
     u'/u/xavier/Dropbox//QSOPairs/data/ESI_redux/SDSSJ002802.61-104936.1_F.fits.gz']



Meta
----

.. code:: python

    reload(pbuild)
    meta = pbuild.mk_meta(flux_files, fname=True, skip_badz=True)


.. parsed-literal::

    Using /u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 for the catalog file
    Using /u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 for the DB file
    Available surveys: [u'GGG']


.. parsed-literal::

    /Users/xavier/local/Python/igmspec/igmspec/privatedb/build.py:138: UserWarning: Skipping 12 entries without a parseable redshift
      np.sum(badz)))


Spectra
-------

.. code:: python

    hdf = h5py.File('tmp.hdf5','w')

.. code:: python

    reload(pbuild)
    pbuild.ingest_spectra(hdf, 'test', meta)


.. parsed-literal::

    Adding test survey to DB


.. parsed-literal::

    /Users/xavier/local/Python/igmspec/igmspec/ingest/utils.py:57: UserWarning: unicode in column SPEC_FILE.  Will convert to str for hdf5
      warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))
    /Users/xavier/local/Python/igmspec/igmspec/ingest/utils.py:57: UserWarning: unicode in column GRATING.  Will convert to str for hdf5
      warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))
    /Users/xavier/local/Python/igmspec/igmspec/ingest/utils.py:57: UserWarning: unicode in column INSTR.  Will convert to str for hdf5
      warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))
    /Users/xavier/local/Python/igmspec/igmspec/ingest/utils.py:57: UserWarning: unicode in column TELESCOPE.  Will convert to str for hdf5
      warnings.warn("unicode in column {:s}.  Will convert to str for hdf5".format(key))


.. code:: python

    hdf.close()

All in One
----------

.. code:: python

    from igmspec.privatedb import build as pbuild

.. code:: python

    tree = os.getenv('DROPBOX_DIR')+'/QSOPairs/data/MMT_redux/'

.. code:: python

    reload(pbuild)
    pbuild.mk_db([tree], ['test'], 'tmp.hdf5',fname=True, skip_badz=True)


.. parsed-literal::

    Using /u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 for the catalog file
    Using /u/xavier/local/Python/igmspec/DB/IGMspec_DB_v01.hdf5 for the DB file
    Available surveys: [u'GGG']
    Adding test survey to DB
    Wrote tmp.hdf5 DB file


