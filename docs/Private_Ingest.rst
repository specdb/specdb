
Ingesting Private Datasets (v2.0)
=================================

:download:`Download <nb/Private_Ingest.ipynb>` this notebook.

.. code:: python

    # imports
    import warnings
    warnings.filterwarnings('ignore')
    
    import h5py
    import specdb
    import glob
    
    from astropy.table import Table
    from linetools import utils as ltu
    
    from specdb.build import privatedb as pbuild
    from specdb.build import utils as spbu
    from specdb.specdb import IgmSpec

Test on Single Folder
---------------------

.. code:: python

    tree = specdb.__path__[0]+'/build/tests/files'
    #os.getenv('DROPBOX_DIR')+'/QSOPairs/data/MMT_redux/'

.. code:: python

    reload(pbuild)
    flux_files = pbuild.grab_files(tree)
    len(flux_files)




.. parsed-literal::

    3



.. code:: python

    flux_files[:5]




.. parsed-literal::

    ([u'/Users/xavier/local/Python/specdb/specdb/build/tests/files//./SDSSJ001605.89+005654.3_b800_F.fits.gz',
      u'/Users/xavier/local/Python/specdb/specdb/build/tests/files//./SDSSJ001607.27+005653.1_b800_F.fits.gz'],
     None,
     None)



--------------

Directory Tree -- Step by Step
------------------------------

.. code:: python

    tree2 = specdb.__path__[0]+'/data/test_privateDB'

.. code:: python

    branches = glob.glob(tree2+'/*')
    branches[0]




.. parsed-literal::

    '/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS'



Get started
~~~~~~~~~~~

.. code:: python

    id_key = 'TEST_ID'
    maindb, tkeys = spbu.start_maindb(id_key)

Files
~~~~~

.. code:: python

    reload(pbuild)
    mflux_files, meta_file, _ = pbuild.grab_files(branches[0])
    len(mflux_files)




.. parsed-literal::

    2



.. code:: python

    mflux_files[:5]




.. parsed-literal::

    [u'/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095240.17+515250.03.fits.gz',
     u'/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095243.05+515121.15.fits.gz']



.. code:: python

    meta_file




.. parsed-literal::

    u'/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS/COS_meta.json'



.. code:: python

    meta_dict = ltu.loadjson(meta_file)
    meta_dict




.. parsed-literal::

    {u'maxpix': 60000,
     u'meta_dict': {u'TELESCOPE': u'HST'},
     u'parse_head': {u'DATE-OBS': u'DATE',
      u'GRATING': u'OPT_ELEM',
      u'INSTR': u'INSTRUME',
      u'R': True}}



ztbl (read from file)
~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ztbl = Table.read(specdb.__path__[0]+'/data/test_privateDB/testDB_ztbl.fits')
    ztbl




.. raw:: html

    &lt;Table length=6&gt;
    <table id="table4652212112" class="table-striped table-bordered table-condensed">
    <thead><tr><th>RA</th><th>DEC</th><th>ZEM</th><th>ZEM_SOURCE</th><th>SPEC_FILE</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>str5</th><th>str35</th></tr></thead>
    <tr><td>331.992916667</td><td>12.9956388889</td><td>1.0</td><td>UNKNW</td><td>SDSSJ220758.30+125944.3_F.fits</td></tr>
    <tr><td>261.35275</td><td>30.6344166667</td><td>1.1</td><td>UNKNW</td><td>SDSSJ172524.66+303803.9_F.fits</td></tr>
    <tr><td>345.184833333</td><td>1.92825</td><td>1.2</td><td>UNKNW</td><td>SDSSJ230044.36+015541.7_r600_F.fits</td></tr>
    <tr><td>345.184833333</td><td>1.92825</td><td>1.2</td><td>UNKNW</td><td>SDSSJ230044.36+015541.7_b400_F.fits</td></tr>
    <tr><td>148.167375</td><td>51.8805638889</td><td>1.3</td><td>UNKNW</td><td>J095240.17+515250.03.fits.gz</td></tr>
    <tr><td>148.179375</td><td>51.855875</td><td>1.4</td><td>UNKNW</td><td>J095243.05+515121.15.fits.gz</td></tr>
    </table>



Meta
~~~~

.. code:: python

    reload(pbuild)
    meta = pbuild.mk_meta(mflux_files, ztbl, fname=True, mdict=meta_dict['meta_dict'], parse_head=meta_dict['parse_head'])

.. code:: python

    meta[0:3]




.. raw:: html

    &lt;Table length=2&gt;
    <table id="table4652211664" class="table-striped table-bordered table-condensed">
    <thead><tr><th>RA_GROUP</th><th>DEC_GROUP</th><th>STYPE</th><th>zem_GROUP</th><th>sig_zem</th><th>flag_zem</th><th>SPEC_FILE</th><th>DATE-OBS</th><th>GRATING</th><th>R</th><th>INSTR</th><th>TELESCOPE</th><th>EPOCH</th><th>GROUP_ID</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>str3</th><th>float64</th><th>float64</th><th>str8</th><th>unicode96</th><th>str10</th><th>str5</th><th>float64</th><th>str3</th><th>unicode3</th><th>float64</th><th>int64</th></tr></thead>
    <tr><td>148.167375</td><td>51.8805638889</td><td>QSO</td><td>1.3</td><td>0.0</td><td>UNKNW</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095240.17+515250.03.fits.gz</td><td>2015-05-31</td><td>G130M</td><td>17000.0</td><td>COS</td><td>HST</td><td>2000.0</td><td>0</td></tr>
    <tr><td>148.179375</td><td>51.855875</td><td>QSO</td><td>1.4</td><td>0.0</td><td>UNKNW</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095243.05+515121.15.fits.gz</td><td>2015-12-08</td><td>G130M</td><td>17000.0</td><td>COS</td><td>HST</td><td>2000.0</td><td>1</td></tr>
    </table>



Without fname=True
^^^^^^^^^^^^^^^^^^

::

    Requires SPEC_FILE in ztbl

.. code:: python

    meta2 = pbuild.mk_meta(mflux_files, ztbl, fname=False, mdict=meta_dict['meta_dict'], parse_head=meta_dict['parse_head'])

.. code:: python

    meta2




.. raw:: html

    &lt;Table length=2&gt;
    <table id="table4653773904" class="table-striped table-bordered table-condensed">
    <thead><tr><th>RA_GROUP</th><th>DEC_GROUP</th><th>STYPE</th><th>zem_GROUP</th><th>sig_zem</th><th>flag_zem</th><th>SPEC_FILE</th><th>DATE-OBS</th><th>GRATING</th><th>R</th><th>INSTR</th><th>TELESCOPE</th><th>EPOCH</th><th>GROUP_ID</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>str3</th><th>float64</th><th>float64</th><th>str8</th><th>unicode96</th><th>str10</th><th>str5</th><th>float64</th><th>str3</th><th>unicode3</th><th>float64</th><th>int64</th></tr></thead>
    <tr><td>148.167375</td><td>51.8805638889</td><td>QSO</td><td>1.3</td><td>0.0</td><td>UNKNW</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095240.17+515250.03.fits.gz</td><td>2015-05-31</td><td>G130M</td><td>17000.0</td><td>COS</td><td>HST</td><td>2000.0</td><td>0</td></tr>
    <tr><td>148.179375</td><td>51.855875</td><td>QSO</td><td>1.4</td><td>0.0</td><td>UNKNW</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS//./J095243.05+515121.15.fits.gz</td><td>2015-12-08</td><td>G130M</td><td>17000.0</td><td>COS</td><td>HST</td><td>2000.0</td><td>1</td></tr>
    </table>



Add Group and IDs
~~~~~~~~~~~~~~~~~

.. code:: python

    gdict = {}
    flag_g = spbu.add_to_group_dict('COS', gdict)
    maindb = pbuild.add_ids(maindb, meta, flag_g, tkeys, id_key, first=(flag_g==1))


.. parsed-literal::

    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... INSTR TELESCOPE EPOCH GROUP_ID
    -------- --------- ----- --------- ------- ... ----- --------- ----- --------


.. code:: python

    maindb




.. raw:: html

    &lt;Table length=2&gt;
    <table id="table4585864144" class="table-striped table-bordered table-condensed">
    <thead><tr><th>flag_group</th><th>sig_zem</th><th>flag_zem</th><th>RA</th><th>DEC</th><th>STYPE</th><th>zem</th><th>TEST_ID</th></tr></thead>
    <thead><tr><th>int64</th><th>float64</th><th>str8</th><th>float64</th><th>float64</th><th>str3</th><th>float64</th><th>int64</th></tr></thead>
    <tr><td>1</td><td>0.0</td><td>UNKNW</td><td>148.167375</td><td>51.8805638889</td><td>QSO</td><td>1.3</td><td>0</td></tr>
    <tr><td>1</td><td>0.0</td><td>UNKNW</td><td>148.179375</td><td>51.855875</td><td>QSO</td><td>1.4</td><td>1</td></tr>
    </table>



.. code:: python

    gdict




.. parsed-literal::

    {'COS': 1}



Spectra
~~~~~~~

.. code:: python

    hdf = h5py.File('tmp.hdf5','w')

.. code:: python

    reload(pbuild)
    pbuild.ingest_spectra(hdf, 'test', meta, max_npix=meta_dict['maxpix'])


.. parsed-literal::

    Adding test group to DB


Finish
~~~~~~

.. code:: python

    pbuild.write_hdf(hdf, 'TEST_DB', maindb, [str('SDSS')], gdict, 'v01')

Directory Tree -- All in One
----------------------------

.. code:: python

    ztbl = Table.read(specdb.__path__[0]+'/data/test_privateDB/testDB_ztbl.fits')

.. code:: python

    reload(pbuild)
    pbuild.mk_db('TEST_DB', tree2, 'tmp.hdf5', ztbl, fname=True)


.. parsed-literal::

    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... INSTR TELESCOPE EPOCH GROUP_ID
    -------- --------- ----- --------- ------- ... ----- --------- ----- --------
    Adding COS group to DB
    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/ESI
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... GRATING EPOCH GROUP_ID tGRB
    -------- --------- ----- --------- ------- ... ------- ----- -------- ----
    Adding ESI group to DB
    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/LRIS
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... GRATING TELESCOPE EPOCH GROUP_ID
    -------- --------- ----- --------- ------- ... ------- --------- ----- --------
    Adding LRIS group to DB
    Wrote tmp.hdf5 DB file


.. code:: python

    # Without fname
    pbuild.mk_db('TEST_DB', tree2, 'tmp2.hdf5', ztbl, fname=False)


.. parsed-literal::

    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/COS
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... INSTR TELESCOPE EPOCH GROUP_ID
    -------- --------- ----- --------- ------- ... ----- --------- ----- --------
    Adding COS group to DB
    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/ESI
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... GRATING EPOCH GROUP_ID tGRB
    -------- --------- ----- --------- ------- ... ------- ----- -------- ----
    Adding ESI group to DB
    Working on branch: /Users/xavier/local/Python/specdb/specdb/data/test_privateDB/LRIS
    The following sources were previously in the DB
    RA_GROUP DEC_GROUP STYPE zem_GROUP sig_zem ... GRATING TELESCOPE EPOCH GROUP_ID
    -------- --------- ----- --------- ------- ... ------- --------- ----- --------
    Adding LRIS group to DB
    Wrote tmp2.hdf5 DB file


By script
~~~~~~~~~

::

    specdb_privatedb testDB ../../specdb/data/test_privateDB tst3_DB.hdf5

Check ESI meta
~~~~~~~~~~~~~~

.. code:: python

    igmsp = IgmSpec(db_file='tmp2.hdf5', verbose=True)


.. parsed-literal::

    Using tmp2.hdf5 for the DB file
    Available groups: [u'COS', u'ESI', u'LRIS']


.. code:: python

    igmsp['ESI'].meta




.. raw:: html

    &lt;Table length=2&gt;
    <table id="table4666987408" class="table-striped table-bordered table-condensed">
    <thead><tr><th>RA_GROUP</th><th>DEC_GROUP</th><th>STYPE</th><th>zem_GROUP</th><th>sig_zem</th><th>flag_zem</th><th>DATE-OBS</th><th>R</th><th>EPOCH</th><th>GROUP_ID</th><th>tGRB</th><th>PRIV_ID</th><th>NPIX</th><th>WV_MIN</th><th>WV_MAX</th><th>SPEC_FILE</th><th>INSTR</th><th>TELESCOPE</th><th>GRATING</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>str3</th><th>float64</th><th>float64</th><th>str8</th><th>str10</th><th>float64</th><th>float64</th><th>int64</th><th>str21</th><th>int64</th><th>int64</th><th>float64</th><th>float64</th><th>str98</th><th>str3</th><th>str7</th><th>str3</th></tr></thead>
    <tr><td>261.3528</td><td>30.6344</td><td>QSO</td><td>1.100</td><td>0.0</td><td>UNKNW</td><td>2015-05-19</td><td>4545.0</td><td>2000.0</td><td>0</td><td>2009-11-23:10:12:13.2</td><td>2</td><td>27931</td><td>3993.5</td><td>10131.6</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/ESI//./SDSSJ172524.66+303803.9_F.fits</td><td>ESI</td><td>Keck-II</td><td>ECH</td></tr>
    <tr><td>331.9929</td><td>12.9956</td><td>QSO</td><td>1.000</td><td>0.0</td><td>UNKNW</td><td>2008-06-04</td><td>4545.0</td><td>2000.0</td><td>1</td><td>2007-08-13:10:22:23.3</td><td>3</td><td>27926</td><td>3993.5</td><td>10129.9</td><td>/Users/xavier/local/Python/specdb/specdb/data/test_privateDB/ESI//./SDSSJ220758.30+125944.3_F.fits</td><td>ESI</td><td>Keck-II</td><td>ECH</td></tr>
    </table>



--------------

JSON files for meta table
-------------------------

.. code:: python

    parse_head = {'DATE-OBS':'DATE', 'TELESCOPE':'TELESCOP','INSTR':'INSTRUME', 'R': True}
    mdict = dict(GRATING='ALL', R=8000.)

.. code:: python

    db_dict = dict(parse_head=parse_head, meta_dict=mdict, maxpix=60000)

.. code:: python

    jdict = ltu.jsonify(db_dict)

.. code:: python

    ltu.savejson('tst.json', jdict, easy_to_read=True, overwrite=True)
