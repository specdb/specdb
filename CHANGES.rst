specdb Changelog
================

Unreleased (2026)
-----------------

This release modernizes the package after four dormant years.  Before it, the
test suite did not run to completion on any current scientific-Python stack.

Breaking changes
^^^^^^^^^^^^^^^^

- **Python 3.11 is now the minimum** (``requires-python = ">=3.11"``).  Users on
  older Pythons will silently receive an older release from PyPI rather than an
  error.
- **``setup.py`` and ``setup.cfg`` are gone**, replaced by a PEP 621
  ``pyproject.toml``.  ``python setup.py install`` no longer works; use
  ``pip install .``.
- **The ``bin/`` scripts are gone**, replaced by ``console_scripts`` entry
  points of the same names (``specdb_cat``, ``specdb_chk``,
  ``specdb_get_igmspec``, ``specdb_get_qpq``, ``specdb_get_uvqs``,
  ``specdb_meta``, ``specdb_plot``, ``specdb_privatedb``, ``specdb_sdss``).
  Each ``specdb.scripts.*.main()`` now takes ``args=None`` and calls
  ``parser()`` itself when called with no arguments.
- **Debugger calls have been replaced by exceptions.**  39 ``pdb.set_trace()``
  calls in shipped modules now raise ``ValueError``/``IndexError``/
  ``NotImplementedError`` with a message.  Code that previously dropped the
  user into an interactive debugger (or hung, under pytest or in CI) now fails
  loudly.  Affected in particular: ``QueryCatalog.pairs()`` with ``dv <= 0``
  raises ``NotImplementedError``; ``mk_db`` refuses to write a database whose
  meta data will not stack; ``ingest_spectra(chk_meta_only=True)`` now returns
  instead of stopping.
- ``mk_meta``'s ``chkz`` argument is deprecated.  It existed only to open a
  debugger; a ``ValueError`` is now always raised for unparseable redshifts.
- ``zem_from_radec``'s ``debug`` argument no longer opens a debugger.
- ``specdb.setup_package.py`` (an unused astropy-helpers leftover) was removed.

Bug fixes
^^^^^^^^^

- **``np.in1d`` was removed in numpy 2.0.**  Replaced with ``np.isin`` at eight
  sites.  This alone had broken ``query_dict``, ``chk_in_group``,
  ``find_ids_in_groups``, ``spectra_from_*``, ``meta_from_*``, ``get_sdss``,
  ``mk_db``, ``ingest_spectra`` and every command-line script.
- **``astropy.table.operations.common_dtype`` moved** to
  ``astropy.utils.metadata`` in astropy 8, breaking ``clean_vstack``.  A bare
  ``except:`` in ``build.utils.chk_vstack`` had been swallowing the error, so
  every private-database build reported "Meta data will not stack" with no
  diagnostic.
- **SSA/VOTable generation was entirely broken under astropy 8.**  The
  hand-rolled ``{'version_1_3_or_later': True}`` config lacked the
  ``version_1_4_or_later`` key astropy now requires; the config is now taken
  from the ``VOTableFile`` itself.
- ``QueryCatalog.cat_from_coords`` called ``self.match_coord``, which was
  renamed to ``query_coords`` in 2017 and left this one caller broken.
- ``query_coords`` raised ``TypeError`` for a length-1 ``SkyCoord`` *array*
  (``int()`` on a one-element ndarray is an error in numpy 2), and built a
  ``(1, 1)`` ID array.  Scalar, length-1 and N-element input now all behave.
- ``cat_utils.write_cat_to_fits`` raised ``ValueError``: ``Table()`` will not
  accept a bare ``h5py.Dataset``.
- ``build.utils.parse_header`` corrupted T-separated ISO timestamps: it
  replaced every ``:`` with ``-``, turning ``2015-09-10T08:36:38`` into
  ``2015-09-10T08-36-38``.  astropy now sniffs the format first.
- ``build.privatedb.grab_files`` returned ``glob.glob`` output unsorted, so
  meta-data row order depended on the filesystem.  Now sorted.
- ``zem.utils.zem_from_radec`` built its ``zsource`` array as ``<U8``, silently
  truncating any longer name.  ``BOSS_PIPE`` -- a real ``z_priority`` entry --
  became ``BOSS_PIP`` in **every private database built with ``mk_meta``**.
- ``utils.hdf_decode`` raised ``RuntimeError`` on a bytes-keyed dict (it
  mutated the dict while iterating it) and left the undecoded keys behind.
- ``SpecDB.__init__`` reported "DB not found" for any failure at all, including
  a broken ``igmspec`` install.  The original exception is now chained.
- The COS branch of ``set_resolution`` reported the wrong header card
  (``DECKNAME``) when ``OPT_ELEM`` was missing.
- ``specdb.py`` used ``verbose & (...)``, which prints nothing for any truthy
  ``verbose`` other than ``True``.  Now ``and``.
- Integer type checks now use ``np.integer`` rather than ``int``/``np.int64``,
  which missed ``np.int32`` and ``np.intp`` on other platforms.

Testing and infrastructure
^^^^^^^^^^^^^^^^^^^^^^^^^^

- **``specdb/build/tests/`` was never collected.**  pytest's default
  ``norecursedirs`` contains ``build``, which matched the ``specdb/build``
  subpackage; seven tests had been invisible since the directory was created.
  ``norecursedirs`` is now set explicitly.  Collection went from 35 to 42.
- The test suite no longer writes into the source tree or the installed
  package.  ``tmp.hdf5``, ``tst_db.hdf5`` and ``tst.xml`` now go to
  ``tmp_path``.  These had been shipping inside built wheels.
- Three assertions that could not fail were repaired
  (``assert IDs[2] == IDs[2]``, a self-comparison of ``GROUP_ID``, and a
  comparison with no ``assert`` at all).  All of the repaired assertions pass,
  so the code under them was correct; the assertions were merely worthless.
- ``test_clean_vstack``'s ``IGMSPEC_DB`` gate now reads ``SPECDB``, the variable
  the library actually uses.
- Added tests for ``cat_from_coords``, single-element-array ``query_coords``
  and ``write_cat_to_fits``.
- Removed a test body disabled by wrapping it in a string literal, which
  referenced the long-removed ``coords_to_spectra``.
- **Continuous integration exists again**: a GitHub Actions workflow covering
  Python 3.11-3.14 on Linux plus macOS, and a ``sphinx-build -W`` docs job.
  ``.travis.yml`` (which cloned the archived ``astropy/ci-helpers`` over the
  dead ``git://`` protocol and ran the removed ``python setup.py test``) was
  deleted.

Documentation
^^^^^^^^^^^^^

- ``sphinx-build`` failed outright under Sphinx 8+ on a legacy
  ``intersphinx_mapping`` form.  The docs now build clean with ``-W``
  (previously 42 warnings, then 7, now 0).
- Added ``.readthedocs.yaml``; without it RTD cannot rebuild, which is why the
  published documentation is still a 2022 Sphinx 1.8.6 build.
- 26 ``.. highlight::`` directives named pages rather than Pygments lexers.
- ``docs/installing.rst`` advertised Python 2.7 and numpy 1.12.
- ``docs/nb_to_rst.src`` had a broken shebang (``#/bin/bash``), a macOS-only
  ``sed -i ''`` and converted a notebook that does not exist.

Python-2 scaffolding removed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

25 ``from __future__`` imports, 10 ``# TEST_UNICODE_LITERALS`` markers, three
``basestring`` shims, four unused ``ustr`` shims, two ``bstr`` shims and six
``str('literal')`` wrappings.  The ``IPython`` ``embed`` imports became unused
once the debugger calls were removed, so IPython is no longer a hard runtime
dependency of ``QueryCatalog``.

Known open questions
^^^^^^^^^^^^^^^^^^^^

The licence is unresolved: ``LICENSE`` is the verbatim GPLv3 with an unfilled
copyright line, while the old ``setup.py`` declared BSD.  ``pyproject.toml``
therefore declares no licence at all.  See ``claude_prompts/2026_update_prompt.md``.
