# 2026 Update

## Goals

This document modernizes `specdb`: it repairs a test suite that **does not
currently run to completion**, replaces packaging that predates PEP 621,
restores continuous integration (there is none), and removes a large amount of
Python-2-era and debugger scaffolding.

`specdb` is a public, multi-author package with downstream users (`igmspec`,
`qpq`, `uvqs`, `casbah`, and anything built on a private specdb database).
Preserve the public API, deprecate rather than delete, and note any behaviour
change so it can go into `CHANGES.rst` -- **which does not exist yet; prompt #7
creates it.**

Run Python only via `conda run -n astro`.  **Xavier performs all git commands.**

Each numbered prompt below is meant to be executed in its own session.  Do only
the task you are pointed at, then append a dated entry under `## Logs`.

## Prompts

1. Read this file.  Execute the 1st task under "Unbreak the suite: `np.in1d`"
2. Read this file.  Execute the 1st task under "Debugger calls in library code"
3. Read this file.  Execute the 1st task under "The uncollected test directory"
4. Read this file.  Execute the 1st task under "astropy 8 compatibility"
5. Read this file.  Execute the 1st task under "SSA and VOTable"
6. Read this file.  Execute the 1st task under "The license contradiction"
7. Read this file.  Execute the 1st task under "Migrate to pyproject.toml"
8. Read this file.  Execute the 1st task under "Build CI from scratch"
9. Read this file.  Execute the 1st task under "README, docs and Read the Docs"
10. Read this file.  Execute the 1st task under "Test hygiene"
11. Read this file.  Execute the 1st task under "Python-2 scaffolding"
12. Read this file.  Execute the 1st task under "Bytes and str in the catalog"

**Order matters far more here than it did for `linetools`.**  Prompts 1, 2 and
3 must run first and in that order: until they are done there is no usable
baseline, and no later prompt can tell whether it broke something.  7 -> 8
are sequenced (CI depends on the packaging metadata).  6 blocks any release
but nothing else.  4, 5, 9, 10, 11, 12 are independent of one another once
1-3 have landed.

13. Read this file.  Check the answers I have given under "Q&A".  Proceed
    accordingly and report back.

### Baseline -- the suite is RED

Unlike `linetools`, **there is no green baseline to protect.**  Establish it
like this and read the numbers carefully:

```
conda run -n astro pytest specdb -q
```

As of 2026-09-27 this reports **8 failed, 6 passed** and then *stops*:

```
>>>>>>>>>>>>>>>>>>> PDB set_trace (IO-capturing turned off) >>>>>>>>>>>>>>>>>>>>
> /Users/xavier/Astronomy/specdb/specdb/interface_group.py(183)ids_to_allrows()
-> pdb.set_trace()
!!!!!!!!!!!!!!!!!!! _pytest.outcomes.Exit: Quitting debugger !!!!!!!!!!!!!!!!!!!
```

Only **14 of 35** collected tests ever execute.  The run does not fail, it
*aborts* -- library code drops into an interactive debugger and pytest quits.
This is also what would happen in CI: a hang or an abort, not a test report.

To see the whole picture before prompt #2 lands, disable the debugger plugin:

```
conda run -n astro pytest specdb -q -p no:debugging < /dev/null
conda run -n astro pytest specdb/build/tests -q -p no:debugging < /dev/null
```

which gives **23 failed / 12 passed** and **5 failed / 2 passed** respectively.

**The true baseline is therefore 42 tests: 28 failed, 14 passed, 0 skipped.**
(`specdb/build/tests` is not collected by `pytest specdb` at all -- see
prompt #3.)  Every later prompt must state the new counts explicitly and
account for any difference.

The `astro` environment as measured on 2026-09-27: Python 3.14.6, numpy 2.5.2,
astropy 8.0.1, h5py 3.16.0, scipy 1.18.0, matplotlib 3.11.1, pytest 9.1.1,
`linetools` 0.3.3.dev23 installed editable from `/Users/xavier/Astronomy/linetools`.
`specdb` itself is **not** pip-installed; it imports from the checkout.  The
last commit on `master` is `def79be`, 2022-05-23.

---

## Unbreak the suite: `np.in1d`

1. **`np.in1d` was removed in numpy 2.0.**  It is the single cause of 23 of the
   28 failures.  `specdb` calls it at eight sites:

   | Site | Context |
   |---|---|
   | `specdb/cat_utils.py:33` | `match_ids()` -- the shared ID matcher |
   | `specdb/interface_group.py:142` | `ids_to_firstrow()` |
   | `specdb/interface_group.py:146` | `ids_to_firstrow()` |
   | `specdb/interface_group.py:181` | `ids_to_allrows()` -- inside a bare `except` |
   | `specdb/interface_group.py:187` | `ids_to_allrows()` |
   | `specdb/query_catalog.py:536` | `cat_from_ids()` |
   | `specdb/utils.py:219` | `query_table()`, the list/str branch |
   | `specdb/build/utils.py:231` | `chk_meta()` |

   Replace each with `np.isin`.  I checked the inputs at all eight sites: they
   are 1-D ID arrays or astropy `Column.data`, so `np.isin` is a true drop-in
   there.  **Know the difference anyway** -- `np.in1d` always returns a flat
   array of `len(ar1.ravel())`, whereas `np.isin` returns something with the
   shape of `ar1`.  If you find a site where `ar1` is not 1-D, stop and say so
   rather than silently changing the returned shape.

2. This is the fix that makes almost the entire public API work again.
   `qcat.query_dict`, `chk_in_group`, `find_ids_in_groups`, `spectra_from_*`,
   `meta_from_*`, `get_sdss`, `mk_db`, `ingest_spectra` and every CLI script
   under `bin/` currently raise `AttributeError` on any modern numpy.

3. Do **only** this in prompt #1.  It is a mechanical eight-line change and it
   should land on its own so the diff is reviewable and the improvement in the
   test counts is unambiguous.  Expect the failure count to drop from 28 to
   roughly 5; state the actual numbers.

4. Note in `CHANGES.rst` (prompt #7 creates the file; if you are running this
   prompt first, keep a note for whoever does) that `specdb` now requires
   numpy >= 1.13 for `np.isin`, which is not a real constraint.

---

## Debugger calls in library code

1. **There are 39 `pdb.set_trace()` calls in shipped `specdb` modules**, plus 2
   more in `specdb/tests/mk_debug_files.py` (41 total), and `import pdb` at
   module scope in 21 files.  Distribution:

   | File | Count |
   |---|---|
   | `specdb/build/utils.py` | 21 |
   | `specdb/build/privatedb.py` | 12 |
   | `specdb/zem/utils.py` | 2 |
   | `specdb/defs.py` | 1 (line 355) |
   | `specdb/interface_group.py` | 1 (line 183) |
   | `specdb/query_catalog.py` | 1 (line 287) |
   | `specdb/utils.py` | 1 (line 214) |
   | `specdb/tests/mk_debug_files.py` | 2 |

   These are not harmless.  A user who calls `QueryCatalog.pairs()` with
   `dv <= 0` gets an interactive debugger instead of an exception.  Under
   pytest they abort the entire session.

2. The worst instance, and the one that breaks the baseline, is
   `specdb/interface_group.py:180-185`:

   ```python
   try:
       chk_survey = np.in1d(IDs, self.meta[self.idkey])
   except:
       pdb.set_trace()
   if np.sum(chk_survey) != IDs.size:
   ```

   A bare `except:` swallowing everything, then a debugger, then a variable
   that is unbound if the `try` failed.  With `pdb` neutralised this becomes
   `UnboundLocalError: cannot access local variable 'chk_survey'`.  Prompt #1
   removes the proximate cause; this prompt removes the pattern.

3. **The conversion rule.**  A `pdb.set_trace()` sitting where an error should
   be reported must become a raised exception with a message, not a silent
   `pass`.  Specifically:

    - `except: pdb.set_trace()` -> let the exception propagate, or catch the
      *specific* exception and `raise` a `ValueError`/`OSError` that says what
      went wrong.  There are **9 bare `except:` clauses**:
      `specdb/interface_group.py:182`, `specdb/specdb.py:43`,
      `specdb/build/privatedb.py:249` (inside a commented block, dead),
      `specdb/build/privatedb.py:359`, `specdb/build/utils.py:222`,
      `specdb/build/utils.py:269`, and
      `specdb/scripts/specdb_chk.py:52,61,94`.
    - `if <bad condition>: pdb.set_trace()` -> `raise ValueError(...)`.
    - A `set_trace()` that is genuinely a developer breadcrumb on a path that
      cannot be reached -> delete it.

   **Do not convert a `set_trace()` into a bare `pass` or a `warnings.warn`.**
   Every one of these marks a place the original author knew was wrong; turning
   it into a warning converts a loud failure into a silent wrong answer.

4. `specdb/specdb.py:43` deserves individual attention.  Its bare `except:`
   catches `ImportError` from the `import igmspec` fallback at line 567 and
   reports `OSError: DB not found` regardless of the real cause.  A user with a
   broken `igmspec` install is told their database is missing.

5. `from IPython import embed` is imported at **module scope** in
   `specdb/query_catalog.py:19`, `specdb/build/privatedb.py:25` and
   `specdb/scripts/specdb_chk.py:6`, and called at
   `specdb/build/privatedb.py:360`.  That makes IPython a hard runtime
   dependency of the core `QueryCatalog` class.  See Q&A #5 before removing it.

6. Removing `import pdb` from the 21 module headers is the easy part; do it
   last, after every `set_trace()` is gone, so the import error tells you if
   you missed one.

7. Expect the suite to complete without aborting after this prompt, and
   `test_getrows_one_source` / `test_getrows_sources` to stop raising
   `OSError: pytest: reading from stdin while output is captured!`.

---

## The uncollected test directory

1. **`specdb/build/tests/` has never been collected by `pytest specdb`.**  The
   seven tests in it have been invisible to every `pytest` and every
   `python setup.py test` invocation since the directory was created.

   The cause is pytest's built-in default `norecursedirs`, which contains
   `build`.  `specdb/build` is a real source subpackage, but pytest does not
   know that and skips it.  Proof:

   ```
   conda run -n astro pytest --collect-only -q specdb | tail -1
   # 35 tests

   conda run -n astro pytest --collect-only -q -o 'norecursedirs=.* *.egg dist' specdb | tail -1
   # 42 tests
   ```

2. Fix it by setting `norecursedirs` explicitly in the pytest configuration.
   Prompt #7 creates `pyproject.toml`; if you run this prompt first, put it in
   a minimal `[tool.pytest.ini_options]` block there and say so, so prompt #7
   knows not to clobber it:

   ```toml
   [tool.pytest.ini_options]
   norecursedirs = [".*", "*.egg", "dist", "docs/_build"]
   testpaths = ["specdb"]
   ```

   Note the deliberate omission of `build` from that list.  Do **not** also add
   a `testpaths` value that excludes anything -- the `linetools` repo had
   exactly that problem and it hid 200 tests.

3. Verify with `conda run -n astro pytest --collect-only -q | tail -1`: bare
   `pytest` from the repo root must now collect 42, and the seven newly-visible
   tests must appear in the report rather than being silently skipped.

4. **Expect newly-visible failures, and do not paper over them.**  Two of the
   seven fail for reasons unrelated to numpy or astropy:

    - `specdb/build/tests/test_utils.py:36` asserts `'unicode' in dtype.name`.
      On numpy 2 the `.name` of a `<U` dtype is `str192`, not `unicode192`.
      This is a test-side Python-2 assumption; fix the test.
    - `specdb/build/tests/test_privatedb.py:33` compares two redshifts against
      a hard-coded order and gets them swapped
      (`ACTUAL [2.597, 2.395]` vs `DESIRED [2.395, 2.597]`).  The real defect
      is in the library: `grab_files()` (`specdb/build/privatedb.py:68,70`)
      returns `glob.glob` output **unsorted**, so the result depends on
      filesystem order.  Sort it in `grab_files` -- do not sort in the test.
      That makes the function deterministic across machines, which is what the
      test was assuming all along.

---

## astropy 8 compatibility

1. **`astropy.table.operations.common_dtype` no longer exists.**
   `specdb/utils.py:32` does
   `from astropy.table import operations, TableMergeError, vstack` and line 48
   calls `operations.common_dtype([...])`.  In astropy 8.0.1 the function lives
   at `astropy.utils.metadata.common_dtype`.  Verify with

   ```
   conda run -n astro python -c "import astropy.table.operations as o; print('common_dtype' in dir(o))"
   ```

   which prints `False`.

   This breaks `clean_vstack()`, reached from `SpecDB.meta_from_coords`
   (`specdb/specdb.py:204`) whenever matches span more than one group -- i.e.
   the common case for any real database.

2. **It is currently hidden.**  `specdb/build/utils.py:267-270` wraps
   `clean_vstack` in a bare `except: chk = False`, so `chk_vstack` reports a
   plain `False` and `mk_db` (`specdb/build/privatedb.py:660-663`) prints
   "Meta data will not stack" and drops into `pdb.set_trace()` for **every**
   private database build.  Anyone who has tried to build a private specdb in
   the last few years has hit this and been given no diagnostic at all.  Fix
   the import; then make that `except` specific so the next such failure is
   visible.

3. **`int()` on a one-element array is now a `TypeError`.**
   `specdb/query_catalog.py:477-479`:

   ```python
   if len(d2d) == 1:  # Annoying array/scalar bit
       IDs = np.array([self.cat[self.idkey][idx]])
       idx = np.array([int(idx)])
   ```

   With numpy 2.5 this raises
   `TypeError: only 0-dimensional arrays can be converted to Python scalars`.
   Reproduce by passing a length-1 array `SkyCoord(ra=[ra0], dec=[dec0], unit='deg')`
   to `query_coords`; a scalar `SkyCoord` and an N>=2 array both work, so this
   is precisely the case the comment calls "annoying".  Use
   `np.atleast_1d(idx)` (or `idx.reshape(1)`) instead of `int(idx)`.

   Affects `meta_from_coords` and `spectra_in_group`.  **No test covers it** --
   add one.

4. **`Table(h5py.Dataset)` is rejected.**  `specdb/cat_utils.py:118` does
   `Table(hdf['catalog'])` and raises
   `ValueError: Data type <class 'h5py._hl.dataset.Dataset'> not allowed to init Table`.
   It needs `Table(hdf['catalog'][()])`.  `write_cat_to_fits` is a public
   function with no test; add one while you are there.

5. **`QueryCatalog.match_coord` does not exist.**
   `specdb/query_catalog.py:112` (`cat_from_coords`) calls `self.match_coord(...)`
   and raises `AttributeError`.  The method was removed or renamed at some
   point and the caller was not updated; the sibling `radial_search` at line 522
   is a stub that only does `raise DeprecationWarning`.  Work out what
   `cat_from_coords` was supposed to call -- most likely `query_coords` -- and
   say in the log what you concluded and on what evidence.  If you cannot
   establish it with confidence, **stop and ask** rather than guessing at the
   semantics of a public method.

---

## SSA and VOTable

1. **Simple Spectral Access is entirely broken under astropy 8.**
   `specdb/ssa.py:397` builds a config dict by hand:

   ```python
   cdict = {'version_1_3_or_later': True}
   ```

   and passes it as `config=` to `Param(...)` at line 466.  astropy now looks up
   `config["version_1_4_or_later"]` (`astropy/io/votable/tree.py:250`), which is
   absent, so the whole thing raises `KeyError: 'version_1_4_or_later'`.

   Both entry points fail: `SSAInterface.querydata(...)` via `metaquery_param()`
   (`ssa.py:135`) and `querydata(FORMAT='METADATA')` (`ssa.py:74` -> `:334`).
   Three tests fail on this.

   The right fix is almost certainly to stop hand-rolling the config dict and
   let astropy supply it, rather than adding a second magic key that the next
   astropy release will also rename.  Investigate how astropy expects a caller
   to construct a standalone `Param` in version 8, and say in the log what the
   supported route is.

2. While you are in this file, the VOTable output emits **21 warnings per
   run**, all from `specdb`'s own construction, all indicating output that does
   not conform to the VOTable spec.  Since the point of SSA is interoperability
   with IVOA clients, these are not cosmetic:

   | Warning | Count | Site |
   |---|---|---|
   | `W03` implicit ID generated from name `INPUT:*` | 10 | `ssa.py:355,362,367,371,375` |
   | `W48` unknown attribute `content` on `INFO` | 5 | `ssa.py:70,85,100,113,123,325` |
   | `W08` `'value'` must be str or bytes | 3 | `ssa.py:70` (`value=1.1`), `:85` (`value=SIZE`, a float) |
   | `W47` missing `arraysize` | 3 | `ssa.py:435` (`TargetName`) |

   `W48` is the interesting one: `content=` is simply not a VOTable `INFO`
   attribute, so whatever those five strings were meant to convey is being
   dropped on the floor for every client.  Work out what attribute was
   intended (`value`? a child `DESCRIPTION`?) and fix it, or -- if the
   information is genuinely not expressible -- delete the argument and say so.

3. `specdb/tests/test_ssa.py:84` writes `tst.xml` into the package directory
   via `data_path()`.  It is not currently observed only because the test dies
   two lines earlier.  Prompt #10 covers this; if you fix SSA first, you will
   be the one who makes that artifact appear, so fix the write at the same time.

4. Three public functions in this module have no test at all: `ssa_defs`,
   `meta_to_ssa_vo`, `build_metaquery`.

---

## The license contradiction

**This is a question, not a task.  Do not "fix" it by choosing a license.**

1. The repository contains four mutually inconsistent statements about its own
   license:

   | Source | Claim |
   |---|---|
   | `LICENSE` (35,141 bytes, added 2016-09-19, never modified) | **GNU General Public License v3**, verbatim |
   | `setup.py:29` | `license = 'BSD'` |
   | `setup.py:2`, `specdb/setup_package.py:1` | "3-clause BSD style license - see LICENSE.rst" |
   | `bin/*:3` (all nine) | "See top-level LICENSE file for Copyright information" |
   | Zenodo record 1069833 | `other-open` |

   There is **no `LICENSE.rst`** anywhere in the repository.  The `LICENSE`
   file contains zero occurrences of `specdb`, `Prochaska` or `Xavier` -- no
   copyright line was ever filled in, and it still carries the stock FSF
   "How to Apply These Terms" appendix.

2. This matters beyond tidiness.  GPLv3 and BSD-3 are not interchangeable, and
   `specdb` sits in an ecosystem (`astropy`, `linetools`, `igmspec`) that is
   uniformly BSD-3.  If the GPL file is authoritative, downstream BSD packages
   that import `specdb` have a licensing problem they do not know about.  If
   the metadata is authoritative, the repository has been shipping the wrong
   license text for ten years.

3. **Relicensing is not a judgement call to make on someone else's behalf, and
   `specdb` has multiple authors.**  Answer Q&A #1 before anything in this
   section is touched.  Once answered, the task is: put the correct license
   text at `LICENSE`, fill in the copyright holder and years, make
   `setup.py`/`pyproject.toml` agree, fix the nine `bin/` headers and
   `specdb/setup_package.py`, and correct the Zenodo metadata.

4. Nothing else in this document depends on the answer -- but **do not cut a
   release until it is resolved.**

---

## Migrate to pyproject.toml

1. The current packaging is a 2016 DESI/`desiutil` template with the DESI parts
   commented out.  Replace `setup.py` and `setup.cfg` with a PEP 621
   `pyproject.toml`.  What is there now:

    - **`setup.py:60` `setup_requires=['pytest-runner']` actively breaks
      things.**  Every `python setup.py <anything>` invocation downloads
      `pytest-runner` from PyPI into a `.eggs/` directory (silently -- `.eggs/`
      is gitignored), and then `python setup.py test` crashes:
      `ModuleNotFoundError: No module named 'pkg_resources'`, because
      `pkg_resources` is gone in setuptools 84.  `python setup.py test` is the
      command `.travis.yml` uses.  Delete this.
    - **`setup.py:61` `tests_require=['pytest']`** -- setuptools already emits
      `UserWarning: Unknown distribution option: 'tests_require'`.  Dead.
    - **`setup.cfg` is two lines**, `[aliases] test=pytest`.  The `setup.py test`
      command was removed in setuptools 72.  Dead; delete the file.
    - **There is no `install_requires` at all** (`setup.py:53` is commented
      out).  `pip install specdb` installs zero dependencies today.
    - **`setup.py:52` `requires = ['Python (>2.7.0)']`** is legacy distutils
      `Requires` metadata that pip has never read.  There is no
      `python_requires` anywhere.
    - **`setup.py:8` `from distutils.extension import Extension`** is unused,
      and `distutils` was removed from the stdlib in Python 3.12 -- it imports
      here only because setuptools shims it.  Delete the import.
    - `setup.py:3` `from __future__ import ...`, and the commented-out
      `desiutil` / `use_2to3` / `cmdclass` / `entry_points` lines, are all dead.

2. **Dependencies.**  There is no existing list to translate, so build one from
   what the code actually imports.  I verified these by grep:

   `numpy`, `astropy`, `h5py`, `linetools`, `matplotlib`, `psutil`
   (`specdb/interface_group.py:6`, module level), `IPython` (module level in
   three files -- see Q&A #5).  `igmspec` is imported lazily at
   `specdb/specdb.py:567` and is genuinely optional.

   Note what is **not** imported anywhere despite being advertised: `scipy` and
   `pyyaml` appear in `.travis.yml:26` and `docs/installing.rst`, and `PyQt5`
   in the docs, but no `specdb` module imports any of them.  They are probably
   transitive `linetools` dependencies.  Do not declare a dependency `specdb`
   does not use; say in the log which you dropped and why.

   `specdb/scripts/get_uvqs.py:50-55` and `get_qpq.py:50-55` shell out to the
   **`wget` binary** and raise `RuntimeError` if it is absent.  That is an
   undeclarable system dependency; note it in the README rather than pretending
   `pip` can satisfy it.  (Replacing it with `requests` or `urllib` would be a
   real improvement, but it is a behaviour change -- propose, do not do.)

3. **Version.**  There are three disagreeing versions today:
   `setup.py:34` says `0.1.dev0`, `docs/conf.py:79-81` says `0.1`, and the
   GitHub tag and Zenodo record say `v1.0`.  `specdb/__init__.py` is **empty**,
   so `specdb.__version__` raises `AttributeError`.

   Pick a single source of truth.  I recommend `setuptools_scm` as used in
   `linetools`, plus an `importlib.metadata` export so the value survives
   without `setuptools_scm` installed at runtime:

   ```python
   from importlib.metadata import PackageNotFoundError, version as _version

   try:
       __version__ = _version("specdb")
   except PackageNotFoundError:   # running from a source tree
       __version__ = "unknown"
   ```

   See Q&A #2 for what the next version number should be.

4. **Package data is currently incomplete.**  `setup.py:71-81` walks only
   `specdb/data` (15 files, 3.0 MB).  It does **not** include
   `specdb/tests/files/` -- which holds the tracked 7.2 MB
   `IGMspec_DB_v02_debug.hdf5` that every test depends on -- nor
   `specdb/build/tests/files/`.  So a wheel installs a package whose test suite
   cannot run.  `include_package_data=True` at line 81 is inert: there is no
   `MANIFEST.in`.  Also `setup.py:77` builds paths with `'/'.join(...)`, which
   is POSIX-only.

   Ironically `specdb/setup_package.py` names the exact missing directory
   (`{'specdb.tests': ['files/*']}`) -- but **nothing references that file**
   (grep across `*.py`, `*.cfg`, `*.toml`, `*.rst`, `*.yml` returns zero hits).
   It is an astropy-helpers convention and astropy-helpers is not installed.
   Delete it.

   See Q&A #6 on whether 10 MB of test data belongs in the wheel at all.

5. **Console scripts.**  The nine files in `bin/` are identical 16-line shims
   (`#!/usr/bin/env python`, then `import specdb.scripts.X as X; X.main(X.parser())`).
   Every target module defines exactly one `parser()` and one `main()`, so they
   convert cleanly to `[project.scripts]` entry points -- except that the entry
   point convention is `main()` with no arguments, while these are
   `main(parser())`.  Give each `main` an `args=None` default that calls
   `parser()` when nothing is passed, then map:

   `specdb_cat` -> `search_cat`, `specdb_chk` -> `specdb_chk`,
   `specdb_get_igmspec` -> `get_igmspec`, `specdb_get_qpq` -> `get_qpq`,
   `specdb_get_uvqs` -> `get_uvqs`, `specdb_meta` -> `grab_meta`,
   `specdb_plot` -> `plot_specdb`, `specdb_privatedb` -> `specdb_privatedb`,
   `specdb_sdss` -> `sdss_spec`.

   While you are there: the docstrings in `bin/` are copy-paste junk.
   `bin/specdb_cat:9`, `bin/specdb_meta:9` and `bin/specdb_plot:9` all say
   *"This script runs PYPIT"*; four others say *"grabs the igmspec database
   file from the great Internet"*.

6. **Create `CHANGES.rst`.**  There is none.  Seed it with the work from
   prompts 1-6 and this one.

7. **Verify by building, not by testing.**  A green test run proves nothing
   about packaging here, exactly as it did not for `linetools`:

   ```
   conda run -n astro python -m build --sdist --wheel --outdir /tmp/specdb_dist
   conda run -n astro python -c "import zipfile; z=zipfile.ZipFile('/tmp/specdb_dist/....whl'); print([n for n in z.namelist() if 'files/' in n or 'data/' in n])"
   ```

   If `build`, `setuptools` or `setuptools_scm` is missing from `astro`, **say
   so and stop** rather than installing without asking.  Do not declare this
   prompt done on the strength of a passing test suite.

8. Knock-on documentation fixes: `docs/installing.rst:74` tells users to run
   `python setup.py install  # or use develop`.  That will no longer exist.

---

## Build CI from scratch

1. **There is no CI.**  There is no `.github/` directory on this branch or on
   `origin/master` -- verified with `git ls-tree -r origin/master`.  The only
   CI artifact is `.travis.yml`, last touched 2017-11-10, and it cannot run:

    - It clones `git://github.com/astropy/ci-helpers.git` (line 119).  The
      `git://` protocol was disabled by GitHub in 2022, and `astropy/ci-helpers`
      is **archived**.
    - `sudo: false` (line 12) is a removed Travis key.
    - It runs `python setup.py test` (lines 24-25), a command setuptools removed.
    - It targets Python 2.7 and 3.6 (lines 21, 52, 82, 92).
    - Its pycodestyle job (lines 109, 115) lints a directory literally named
      `packagename` -- a cookiecutter placeholder nobody ever substituted.
    - Line 127 references `linetools/tests/coveragerc`, paste-over from another
      repository.

   Delete it and write `.github/workflows/ci_tests.yml`.

2. What the workflow needs, derived from what the tests actually require:

    - **Triggers**: push to `master` and `'2026'`, plus `pull_request`.
      Quote `'2026'` -- bare `- 2026` parses as the YAML *integer* 2026 and will
      never match a branch name.  (This exact trap bit the `linetools` work.)
    - **Matrix**: Python 3.11-3.14 on `ubuntu-latest`, plus one
      `macos-latest` job -- Xavier develops on macOS and nothing has ever
      tested it.  See Q&A #3 for the floor.
    - `actions/checkout@v4` with `fetch-depth: 0` (setuptools_scm needs tags)
      and `actions/setup-python@v5`.
    - **`MPLBACKEND=Agg`.**  `.travis.yml:40` set `SETUP_XVFB=True`; no test
      needs a display, and the backend variable is what actually matters.
    - **No network is required.**  I verified this: no test imports `urllib`,
      `requests` or shells out to `wget`.  Every test that needs a database
      uses the committed 19-source fixture
      `specdb/tests/files/IGMspec_DB_v02_debug.hdf5`, and the build tests use
      committed FITS fixtures.  All 28 current failures are API drift, not
      missing data.  CI can therefore run the full suite offline.
    - Run `pytest` directly.  Make sure it collects 42 (prompt #3).
    - A **docs job**: `sphinx-build -W -b html docs docs/_build/html`.  There is
      no CI on the documentation at all, which is how the Sphinx breakage in
      prompt #9 shipped unnoticed.

3. **`linetools` is the awkward dependency.**  `.travis.yml:121` installed it
   from git HEAD (`git+https://github.com/linetools/linetools.git`), i.e. CI
   always tested against unreleased `linetools`.  The `astro` environment here
   likewise has an editable checkout, not a release.  See Q&A #4 -- this is a
   real decision, not a detail, because `linetools` is itself being modernized
   in parallel and a released version compatible with Python 3.11+ may not
   exist yet.

4. Environment variables the package reads, in case a job ever wants them:
   **`SPECDB`** (`specdb/specdb.py:565,634,687,741`) points at the directory
   holding the real `*_DB_*.hdf5` files; when unset, `IgmSpec()` falls back to
   `import igmspec`.  **`IGMSPEC_DB`** is read by exactly one place --
   `specdb/tests/test_utils.py:37` -- and by **nothing in the library**.  That
   is a bug in the test; see prompt #10.

5. Do not add a coverage badge unless you also wire up the upload.  A badge
   that reports stale data is worse than no badge.

---

## README, docs and Read the Docs

1. **Read the Docs is frozen at 2022.**  `specdb.readthedocs.io` is live and
   serving a build of commit `def79be` dated 2022-05-23, generated with
   **Sphinx 1.8.6**.  There is no `.readthedocs.yaml`, which RTD now requires,
   so the project cannot rebuild even if it were triggered.

2. **A local `sphinx-build` fails immediately** under the Sphinx 9.1.0 in
   `astro`:

   ```
   sphinx.errors.ConfigError: Invalid intersphinx_mapping configuration (1 error)
   ```

   from `docs/conf.py:312`: `intersphinx_mapping = {'https://docs.python.org/': None}`,
   the legacy unnamed form, a hard error since Sphinx 8.  Use the named form.
   With that and `language = None` -> `'en'` (`conf.py:88`; `None` is invalid
   for Sphinx >= 5) the build succeeds with 42 warnings.

3. Those 42 warnings, worth fixing while you are in there:

    - **26 x `Pygments lexer name 'X' is not known`**, where X is `catalog`,
      `groups`, `igmspec`, `meta`, `qpq`, `uvqs`.  Someone used
      `.. highlight:: <page name>` as if it were a section marker:
      `catalog.rst:162`, `groups.rst:39,54,90,107`, `igmspec.rst:37,49,57`,
      `meta.rst` (11 sites, lines 58-217), `qpq.rst:46`,
      `query_dict.rst:45-79` (4 sites), `uvqs.rst:30`.
    - `installing.rst:87` references an undefined label `download-scipts` --
      a typo for `download-scripts`, and that label does not exist either.
    - `docs/_static` and `docs/_templates` are referenced by `conf.py:162,56`
      and do not exist.

4. **The API is undocumented.**  `sphinx.ext.autodoc` is enabled
   (`conf.py:33`) but there is not a single `automodule`/`autoclass`/
   `autofunction` directive in any `.rst` file (grep returns zero).  The
   `sys.path.insert(0, os.path.abspath('../specdb'))` at `conf.py:22` is
   therefore inert -- and wrong if it were ever used: it inserts the *package*
   directory, so `import specdb` would resolve to the `specdb/specdb.py`
   module rather than the package.  Adding an API reference is a larger job;
   propose it, do not do it inside this prompt.

5. **`docs/installing.rst` is stale throughout**: line 31 claims support for
   "2.7, or 3.3 or later"; lines 32-37 pin numpy 1.12 / astropy 1.3 /
   scipy 0.18 / matplotlib 1.4 / PyQt5 / h5py 2.6 (installed: 3.16); line 28
   links to continuum.io; lines 15-16 say `linetools` is "under review for
   astropy affiliated status"; line 61's `linetools` install link returns
   **HTTP 404** (it has a trailing slash after `.html`); line 96 cites
   `paper_coming_soon <https://www.arxiv/>`.

6. **`README.md`.**  The markdown is valid (unlike `linetools`, whose `.md`
   file contained reStructuredText).  The Zenodo badge and DOI both resolve
   (record 1069833, "First Release", v1.0, 2017-12-03).  Fix the typo
   `Virutal` -> `Virtual` on line 12, and add a CI badge once prompt #8 exists.
   Note that Zenodo's license field says `other-open`, which is a third
   answer to the question in prompt #6.

7. **`docs/nb_to_rst.src`** is a 25-line helper that converts three notebooks
   to rst.  Its shebang is `#/bin/bash` -- missing the `!`, so it only works
   when invoked explicitly through bash.  It requires macOS `sed -i ''`, and it
   converts `Simple_Usage.ipynb`, which **does not exist** (the seven notebooks
   in `docs/nb/` do not include it).  Either repair it or delete it; it is
   referenced by neither the Makefile nor `conf.py`.

---

## Test hygiene

1. **Tests write into the source tree and into the package directory.**
   Inventory, verified by grep and by running the suite:

   | Site | Writes | Class |
   |---|---|---|
   | `specdb/build/tests/test_privatedb.py:50` | `tmp.hdf5` (82 KB) in the CWD | repo root |
   | `specdb/build/tests/test_privatedb.py:66` | `specdb/build/tests/files/tst_db.hdf5` (801 KB) | **package dir** |
   | `specdb/tests/test_ssa.py:84` | `specdb/tests/files/tst.xml` | **package dir** |
   | `specdb/tests/mk_debug_files.py:32-33` | overwrites the committed fixture `specdb/tests/files/IGMspec_DB_v02_debug.hdf5` | **package dir**, `__main__` only |

   The `data_path()` writes are the worse category: that directory is read-only
   in a normal installation and shared between parallel test runs.  The
   `tst.xml` write is not observed today only because the test dies earlier
   (prompt #5); fixing SSA will make it appear.

   Verified empirically -- `git status --porcelain` is clean before a run and
   afterwards shows `?? specdb/build/tests/files/tst_db.hdf5` and `?? tmp.hdf5`.
   Nothing in `.gitignore` masks them.  **There is an 801 KB
   `specdb/build/tests/files/tst_db.hdf5` sitting untracked in the tree right
   now**, left by the audit runs that produced this document; delete it before
   you start so your own before/after comparison is meaningful.

2. **Nothing in this repository uses `tmp_path` or `tmpdir`.**  Establish the
   pattern cleanly:

   ```python
   def test_something(tmp_path):
       outfil = str(tmp_path / 'tst_db.hdf5')
   ```

   Afterwards, `git status` must be clean after a full run.  Verify that
   explicitly and say so in the log.  **Do not extend `.gitignore` instead** --
   that is exactly how these went unnoticed.

3. **Dead assertions.**  Three sites assert nothing:

    - `specdb/build/tests/test_utils.py:53`: `assert IDs[2] == IDs[2]` is a
      tautology.  Rows 2 and 3 of that fixture are duplicate coordinates
      `(6., 7.)`, so the intent was almost certainly `IDs[2] == IDs[3]`.
      **Make the change and find out.**  If it then fails, report what
      `get_new_ids` actually returns and stop -- do not loosen it back.
    - `specdb/tests/test_meta_query.py:100`:
      `assert np.sum(meta5['GROUP_ID'] == meta5['GROUP_ID']) == 2` is a
      tautology given the `assert len(meta5) == 2` on the line before.  I
      verified it returns 2 for `[0,0]`, `[0,1]`, `[-1,999]` and `[7,7]` alike.
      Work out which two things were meant to be compared.
    - `specdb/tests/test_meta_query.py:68`: `meta['GROUP'] == 'BOSS_DR12'` with
      **no `assert`** at all -- the result is computed and discarded.

4. Four tests contain no assertion of any kind and only check that nothing
   raises: `test_group_utils.py::test_show_meta`,
   `test_utils.py::test_clean_vstack`, `test_scripts.py::test_plot_spec`,
   `test_scripts.py::test_sdss`.  That is weak but not wrong; strengthen them
   if it is cheap, and note which you left alone.

5. **`specdb/tests/test_utils.py:37` checks the wrong environment variable.**
   It gates a test on `IGMSPEC_DB`, which **nothing in the library ever reads**;
   the library reads `SPECDB`.  So that branch is dead on every machine.  And
   even with the variable set it would fail: I ran the body, and
   `IgmSpec()` warns "SPECDB not set", tries `import igmspec`, gets
   `ModuleNotFoundError`, and raises `OSError: DB not found` from
   `specdb/specdb.py:44`.  Either gate it on `SPECDB` and make it genuinely
   optional, or delete it -- but say which and why.

6. **`specdb/tests/test_get_spectra.py:102-107` is a test body disabled by
   wrapping it in a `'''...'''` string literal.**  It calls
   `igmsp.coords_to_spectra(coords, 'GGG', all_spec=True)`, and
   `coords_to_spectra` **does not exist anywhere in `specdb`** (grep returns
   nothing).  It is a dead reference to a removed API.  Delete the block.
   There is a second such block at `specdb/tests/mk_debug_files.py:57-64`.

7. **Untested public surface**, verified by AST scan -- record this in the log
   even if you do not fix it, because it is the reason so much of the API drift
   above went undetected:

   | Module | Public names | Referenced by any test |
   |---|---|---|
   | `specdb/defs.py` | 9 | **0** |
   | `specdb/zem/utils.py` | 1 (`zem_from_radec`) | **0** |
   | `specdb/query_catalog.py` | `QueryCatalog` | 0 by name (only via `igmsp.qcat`) |
   | `specdb/cat_utils.py` | `flag_to_groups`, `write_cat_to_fits` | 0 |
   | `specdb/utils.py` | `hdf_decode`, `load_db`, `query_table` | 0 by name |
   | `specdb/ssa.py` | `ssa_defs`, `meta_to_ssa_vo`, `build_metaquery` | 0 |
   | `specdb/build/utils.py` | 10 of 16 | 0 |
   | `specdb/specdb.py` | `UVQS`, `QPQ`, `CASBaH` | never instantiated |

8. **`.gitignore`.**  Line 13 reads `./build/` and matches **nothing** -- git
   does not normalize a leading `./`.  Verified with `git check-ignore -v`.
   Be careful fixing it: `specdb/build/` is a real tracked subpackage, so a
   bare `build/` or `build` pattern would ignore `specdb/build/utils.py`.  Only
   the anchored `/build/` is correct.  Also missing and worth adding:
   `.ipynb_checkpoints/` (there are seven notebooks in `docs/nb/`),
   `.DS_Store`, `.pytest_cache/`.  Note that `*.spec` on line 32 (a PyInstaller
   template leftover) would ignore any legitimately named `*.spec` file.

---

## Python-2 scaffolding

*Do this late.  It is a large, purely mechanical diff and it will make every
other prompt's diff unreadable if it lands first.*

1. **`from __future__ import print_function, absolute_import, division, unicode_literals`**
   in 13 files: `specdb/defs.py:3`, `group_utils.py:3`, `interface_group.py:3`,
   `cat_utils.py:3`, `ssa.py:4`, `utils.py:3`, `specdb.py:3`, `zem/utils.py:3`,
   `build/utils.py:3`, `scripts/specdb_privatedb.py:5`, `scripts/get_uvqs.py:5`,
   `scripts/get_igmspec.py:5`, `scripts/get_qpq.py:5`.  Plus more in test files
   and `setup.py:3`.  All no-ops.

2. **Compatibility shims:**

    - `try: basestring / except NameError: basestring = str` --
      `specdb/utils.py:9-12` (used at line 197), `specdb/ssa.py:18-21`
      (**unused**), `specdb/specdb.py:18-21` (used at line 511 in `__getitem__`).
    - `ustr = unicode` / `str` -- `scripts/specdb_privatedb.py:9-12`,
      `get_uvqs.py:10-13`, `get_igmspec.py:10-13`, `get_qpq.py:10-13`.
      **`ustr` is never used anywhere.**  Delete all four.
    - `bstr = bytes` -- `specdb/defs.py:11-14` (unused in that file),
      `specdb/build/utils.py:22-25` (used at lines 230 and 298 -- keep the
      behaviour, just write `bytes`).

3. **`str('...')` wrapping string literals**, an artifact of `unicode_literals`:
   `defs.py:120-121`, `ssa.py:197-198,262-280`, `zem/utils.py:53`,
   `build/privatedb.py:207,483-488,666`, `build/utils.py:422-426`,
   `specdb.py:348`.  Unwrap them.

4. `#reload(spdbu)` commented at `specdb/query_catalog.py:324`.

5. What is **not** here, so you do not go looking: no `print` statements, no
   `except X, e`, no `iteritems`/`has_key`/`xrange`, no `StringIO`, no
   `urllib2`, no old-style classes, no `six`, no `future`.  All 51 `.py` files
   compile cleanly under Python 3.14.  The Python-3 port was done; only the
   scaffolding was left behind.

6. Two near-misses worth fixing while you are in the area, neither of which is
   a Python-2 issue:

    - `specdb/specdb.py:230`: `if self.verbose & ((ii % 100) == 0):`.  With a
      genuine bool this works, but `verbose=2` gives `2 & True == 0`, so it
      never prints.  Use `and`.
    - `specdb/interface_group.py:207` tests `isinstance(x, (int, np.int64))`,
      which misses `np.int32` and `np.intp` on other platforms; lines 115, 139,
      177 and `query_catalog.py:582` accept only Python `int`.  `np.integer` is
      the portable check.

   There are **no** `~`-on-bool bugs in this repository -- I checked all 20 `~`
   sites and every one is on a numpy bool array, which is correct.  (This was a
   real bug in `linetools`; it is not one here.)

---

## Bytes and str in the catalog

*This one is a genuine inconsistency in the public API, not a cleanup.  Read
Q&A #7 before acting.*

1. **Meta tables come back as `str`; the catalog comes back as `bytes`.**

   Group meta tables are read through `spdbu.hdf_decode(..., itype='Table')`
   (`specdb/interface_group.py:70`), which converts every `S` column to `str`.
   The catalog is read with `Table.read(hdf.filename, path='catalog')`
   (`specdb/query_catalog.py:64`) and is **not** decoded, so `cat['STYPE']` is
   `|S3` and `cat['flag_zem']` is `|S10`.

   It mostly works by accident: astropy's `Column.__eq__` transparently encodes
   a `str` comparand, and `query_table` casts the query value with
   `.astype(tbl[key].dtype)` (`specdb/utils.py:217`).  So
   `cat['STYPE'] == 'QSO'` and `query_dict({'STYPE': 'QSO'})` both behave.

   But raw ndarray comparisons silently return all-`False`:
   `np.array([b'QSO']) == 'QSO'` -> `[False, False]`, and
   `np.isin(S_array, ['QSO'])` -> all `False`.  Any user who writes
   `cat['STYPE'].data == 'QSO'`, or calls `np.unique(cat['flag_zem'])` (which
   returns `b'BOSS_PCA '`), or exports the catalog to FITS or CSV, gets bytes
   and no error.

2. **A data bug falls out of this.**  In the debug database `flag_zem` is
   stored as `b'BOSS_PCA '` -- with a trailing space -- so
   `query_dict({'flag_zem': 'BOSS_PCA'})` matches **zero** rows.  numpy `S`
   comparison does not strip trailing whitespace.  Check whether the real
   `IGMspec` and `QPQ` databases have the same padding; if they do, this has
   been silently returning empty results for years.

3. **`specdb/zem/utils.py:53-54` truncates.**

   ```python
   zsource = np.array([str('NONENONE')]*len(ra))   # dtype '<U8'
   ```

   Any source name longer than 8 characters is silently cut.
   `defs.z_priority()` (`specdb/defs.py:135-154`) contains `BOSS_PIPE`, nine
   characters, which becomes `BOSS_PIP`.  This corrupts `flag_zem` in **every
   private database built with `mk_meta`**.  Use an explicit wide dtype or
   `dtype=object`.  This is a real, user-visible data-corruption bug -- give it
   its own `CHANGES.rst` entry.

4. **`hdf_decode` mutates a dict while iterating.**  `specdb/utils.py:92-96`:
   `hdf_decode({b'a': b'b'})` raises
   `RuntimeError: dictionary changed size during iteration`.  Its only caller
   is `specdb/scripts/specdb_chk.py:47`, and h5py 3 always returns `str`
   attribute names, so the branch is dead today -- but it is one h5py change
   away from being live.  Iterate over `list(d.items())`.

5. Two cosmetic-but-confusing side effects to note in the log, not necessarily
   to fix: `hdf_decode` does `remove_column` then re-adds
   (`specdb/utils.py:89-91`), so **all string columns are moved to the end** of
   every meta table; and `specdb/build/privatedb.py:384` fills a missing `R`
   column with the *string* `'DUMMY'`, so `R` becomes a bytes column in the
   meta table instead of a float.

---

## Q&A

1. **The license is contradictory and I will not resolve it for you.**  The
   `LICENSE` file is the verbatim GPLv3 with an unfilled copyright line;
   `setup.py` says `BSD`; the source headers say "3-clause BSD ... see
   LICENSE.rst" and there is no `LICENSE.rst`; Zenodo says `other-open`.
   `specdb` sits in an all-BSD ecosystem (astropy, linetools), so if the GPL
   file is authoritative, downstream BSD packages importing `specdb` have a
   problem they do not know about.  Which is intended, and what should the
   copyright line read (holder and years)?  `specdb` has multiple authors, so
   this may not be yours alone to decide.

   > *Answer:*

2. **Version number.**  Three disagree: `setup.py` says `0.1.dev0`,
   `docs/conf.py` says `0.1`, the GitHub tag and Zenodo say `v1.0`.
   `specdb.__version__` does not exist.  What should the next release be
   called, and do you want `setuptools_scm` to derive it from tags as
   `linetools` now does?

   > *Answer:*

3. **`requires-python` floor.**  `linetools` went to `>=3.11`.  `specdb`
   depends on `linetools`, so it cannot sensibly be lower.  Confirm `>=3.11`,
   and confirm the CI matrix 3.11-3.14 (dropping 3.15 for now, as we did for
   `linetools`)?

   > *Answer:*

4. **How should CI get `linetools`?**  Travis installed it from git HEAD, and
   the `astro` environment has an editable checkout at
   `/Users/xavier/Astronomy/linetools` (0.3.3.dev23), not a release.  Since
   `linetools` is being modernized in parallel, a PyPI release compatible with
   Python 3.11+ may not exist yet.  Pin CI to `linetools` git `master`, or wait
   for a release and depend on that?

   > *Answer:*

5. **`IPython` and `psutil` as hard dependencies.**  `from IPython import embed`
   is at module scope in `specdb/query_catalog.py:19`,
   `specdb/build/privatedb.py:25` and `specdb/scripts/specdb_chk.py:6`, so
   importing the core `QueryCatalog` requires IPython.  `psutil` is imported at
   module scope in `specdb/interface_group.py:6`.  My inclination is to delete
   the `IPython` imports along with the `pdb` ones in prompt #2 (they are
   debugging tools, not features) and declare `psutil` properly.  Agreed, or do
   you use `embed()` deliberately?

   > *Answer:*

6. **Test data in the distribution.**  `specdb/tests/files/IGMspec_DB_v02_debug.hdf5`
   is 7.2 MB and `specdb/data/test_privateDB/` is 3 MB.  The current
   `package_data` ships the 3 MB and omits the 7.2 MB, which means an installed
   `specdb` cannot run its own tests.  Ship both (a ~10 MB wheel), ship
   neither (tests only run from a checkout), or ship only `specdb/data`?

   > *Answer:*

7. **The catalog bytes/str inconsistency (prompt #12).**  Decoding the catalog
   to `str` on read would make it consistent with the meta tables and fix the
   silent all-`False` comparisons -- but it is a **behaviour change** for any
   downstream code that currently compares against `b'...'`.  Fix it, or
   document the inconsistency and leave it?  Related: do the real IGMspec and
   QPQ databases have the trailing-space padding in `flag_zem` that the debug
   database has?

   > *Answer:*

8. **`QueryCatalog.match_coord` (prompt #4, item 5).**  `cat_from_coords` calls
   a method that does not exist, so that public function has been dead for some
   time.  Do you remember what it was meant to call -- `query_coords`? -- or
   should `cat_from_coords` be removed instead?

   > *Answer:*

9. **`specdb/scripts/get_uvqs.py` and `get_qpq.py` shell out to the `wget`
   binary** and raise `RuntimeError` if it is absent.  `get_igmspec.py` has its
   download code commented out and only prints Google Drive links.  Should
   these be rewritten to use `requests`/`urllib` (removing a system dependency
   and making them testable), or left alone?

   > *Answer:*

10. **Scope check.**  This document is considerably larger than the `linetools`
    one because `specdb` has been dormant since 2022-05-23 and the suite does
    not run.  If you would rather stage it -- prompts 1-5 now to get a green
    suite, packaging and CI later -- say so and I will split the file.

    > *Answer:*

### Round 2 -- raised 2026-09-27 while executing prompts 1-12

11. **Prompt #6 (the licence) was not executed.**  Its own text says "This is a
    question, not a task.  Do not 'fix' it by choosing a licence", and Q&A #1 is
    unanswered, so I left `LICENSE` untouched.  One consequence you should know
    about: I deliberately omitted the `license` field from `pyproject.toml`, but
    setuptools globs `LICENSE*` by default, so **the built wheel carries
    `License-File: LICENSE` -- i.e. it ships the GPLv3 text** while declaring no
    licence in its metadata.  That is the least-wrong state I could leave it in,
    but it is not a state to release from.

    > *Answer:*

12. **Prompt #12 item 1 (decoding the catalog) was not executed** either, for
    the same reason -- Q&A #7 is unanswered and it is a behaviour change for
    downstream code comparing against `b'...'`.  I did fix the two unambiguous
    bugs in that section (the `<U8` truncation and the `hdf_decode` dict
    mutation).  I also confirmed the padding: in the debug database `flag_zem`
    really is stored as `b'BOSS_PCA '` with a trailing space, so
    `query_dict({'flag_zem': 'BOSS_PCA'})` matches zero rows.  I cannot check
    the real IGMspec/QPQ databases from here -- they are not on this machine.

    > *Answer:*

13. **Three public keyword arguments now do nothing**, because opening a
    debugger was their entire implementation: `zem_from_radec(debug=...)`,
    `mk_meta(chkz=...)` (which now raises a `DeprecationWarning` if set) and
    the `debug=` breadcrumbs in `build/utils.get_new_ids` / `add_ids`.  I kept
    them in the signatures rather than break the API.  Remove them in a
    deprecation cycle, or give them a real meaning?

    > *Answer:*

14. **`QueryCatalog.pairs(dv <= 0)` now raises `NotImplementedError`.**  The
    `dv > 0` branch (projected pairs) is fully written; the `else` branch was
    only ever `pdb.set_trace()`, so physical pairs were never implemented.  Was
    that deliberate, or is there code elsewhere that should be moved in?

    > *Answer:*

15. **`ingest_spectra(chk_meta_only=True)` now returns instead of stopping.**
    Its docstring says "Only check meta file; will not write", so `return` is
    my reading of the intent -- but the debugger call meant nobody has actually
    exercised this path.  Confirm?

    > *Answer:*

16. **I deleted `bin/`** (all nine shims) when converting them to
    `[project.scripts]` entry points, and gave every `specdb.scripts.*.main()`
    an `args=None` default.  The nine command names are unchanged and I
    verified they install into a clean venv.  Anyone invoking the scripts by
    *path* rather than by name will break.  Fine?

    > *Answer:*

17. **The dependency floors in `pyproject.toml` are my best reading of wheel
    availability for Python 3.11, not a verified fact** (`numpy>=1.24`,
    `astropy>=5.3`, `h5py>=3.8`, `matplotlib>=3.6`, `psutil>=5.9`).  Nothing
    tests them.  `linetools` is declared with no floor at all because it has no
    release that satisfies `requires-python >= 3.11`; CI installs it from git
    master.  This is the same question as Q&A #4 and it now blocks a PyPI
    release: **`pip install specdb` cannot work until `linetools` is released.**

    > *Answer:*

18. **The SSA output still emits four `UnitsWarning`s**: astropy reports that
    `Angstrom` is deprecated in the VOUnit standard and suggests `0.1nm`.
    Changing it is a one-line fix but it changes the VOTable that IVOA clients
    receive, so I left it.  Change it?

    > *Answer:*

19. **`specdb/tests/test_meta_query.py` ends with an unfinished test**: line
    ~127 computes `qdict = dict(DISPERSER='R400')` and then the function ends.
    I left it alone rather than invent an assertion -- I checked, and the
    obvious completion (`meta_from_coords(..., query_dict=qdict)`) is not even a
    valid call, since that method takes `cat_query`/`meta_query`.  What was it
    meant to test?

    > *Answer:*

20. **One correction to this document.**  The "README, docs" section claims the
    `download-scripts` label "does not exist either".  That is wrong: it is at
    `docs/scripts.rst:22`.  The only defect was the `download-scipts` typo in
    the reference, which is fixed.

    > *Answer:*

---

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-27 (Execute prompts 1-12: green suite, PEP 621 packaging, CI, docs)

Executed prompts 1-5 and 7-12.  **Prompt #6 was not executed** -- see Q&A #11.

**Baseline 28 failed / 14 passed (42 tests, run aborting early) -> now 45
passed / 0 failed / 0 skipped.**  The suite grew from 42 to 45 because prompt #4
asked for three new tests.  The docs build clean under `sphinx-build -W`
(previously it did not build at all).  An sdist and a wheel build, and the wheel
installs into a clean venv with all nine console scripts registered.

**1. `np.in1d` -> `np.isin`** at the eight sites listed.  All eight take 1-D ID
arrays or `Column.data`, so the shape difference between the two functions does
not arise.  This one change took the failures from 28 to 9.

**2. Debugger calls.**  39 `pdb.set_trace()` calls removed from shipped modules
(21 in `build/utils.py`, 12 in `build/privatedb.py`, and one each in `defs.py`,
`interface_group.py`, `query_catalog.py`, `utils.py`, plus two in
`zem/utils.py`), 2 more from `tests/mk_debug_files.py`, and `import pdb` from
21 module headers.  Every one became a raised exception with a message, not a
`pass` and not a warning.  The bulk were the `except KeyError: print("Need to
add X"); pdb.set_trace()` pattern in `set_resolution`; those became
`raise ValueError(...)`, which is the error type `set_resolution` already used
elsewhere and the type its caller in `parse_header` already catches.

Nine bare `except:` clauses were narrowed or given a message.  `specdb.py:43`
now chains the original exception (`raise IOError(...) from exc`) instead of
reporting "DB not found" for any failure whatsoever.

Two things I found while doing this that are not in the prompt doc:

- The COS branch of `set_resolution` looked up `head['OPT_ELEM']` but reported
  `head['DECKNAME']` on failure.  Once that became a `raise`, a missing
  `OPT_ELEM` would have raised a *second* `KeyError` inside the handler and
  masked the first.  Fixed to report `OPT_ELEM`.
- **`from IPython import embed` had exactly one call site**, in the bare
  `except:` at `privatedb.py:359`.  Once that was converted, all three
  module-level `embed` imports were dead, so I removed them.  IPython is
  therefore no longer a hard runtime dependency of `QueryCatalog` -- which
  partly answers Q&A #5 without my having to decide anything.

**3. The uncollected directory.**  Confirmed the cause is pytest's built-in
`norecursedirs` containing `build`.  Set it explicitly in `[tool.pytest.ini_options]`
(in a minimal `pyproject.toml`, which prompt #7 then grew).  Collection 35 -> 42.

Of the two newly-visible failures the doc predicted, both were real and both
diagnoses held up.  `test_meta` was **right and the library was wrong**:
`grab_files` returned `glob.glob` output unsorted, so meta-data row order
depended on the filesystem.  Sorting in `grab_files` (both the glob and the
`os.walk` folder list) fixed it, and I did not touch the test.  `test_clean_meta`
was the test's fault -- numpy 2 renamed the `<U` dtype's `.name` from
`unicode192` to `str192` -- and I replaced the string sniffing with
`dtype.kind == 'U'` / `'S'`, which will not drift again.

A **third** newly-visible failure the doc did not anticipate: `test_mkdb` died
in `parse_header` with `Time 2015-09-10T08-36-38 does not match iso format`.
The cause is a genuine bug -- the code replaces *every* `:` with `-` before
parsing, which corrupts the time part of a T-separated ISO timestamp.  It only
shows up for LRIS, whose `parse_head` maps `DATE-OBS` onto the `DATE` card, and
that is the only fixture carrying a full timestamp rather than a bare date.
Now astropy sniffs the format first, with the old munging kept as a fallback.

**4. astropy 8.**  `common_dtype` moved to `astropy.utils.metadata`.  Worth
recording how well hidden this was: `chk_vstack`'s bare `except: chk = False`
turned it into a plain `False`, so `mk_db` printed "Meta data will not stack"
and stopped at a debugger for **every** private-database build, with no
diagnostic at all.  The `except` now prints the exception type and message.

`int()` on a one-element array: the fix in the prompt doc (`np.atleast_1d(idx)`)
was necessary but not sufficient -- with `idx` already a 1-element array,
`np.array([self.cat[self.idkey][idx]])` builds a **(1, 1)** array and the next
line raises `IndexError: too many indices`.  I removed the special case
entirely: `idx` and `d2d` go through `np.atleast_1d`, and `IDs` is built with
`np.array(...)`, which also forces the copy that the `IDs[...] = -1` assignments
below depend on.  Scalar, length-1 and N-element input now all work, and there
is a test for it.

**`match_coord`: the history settles it.**  `git log -S` found commit `a1b4018`
("query list of coords", 2017-01-04), which deletes `def match_coord` and adds
`def query_coords` in the same diff, moving the identical
`match_coordinates_sky(coords, self.coords, nthneighbor=1)` body across.  The
old method returned exactly the ID array with the -1/-2 semantics that
`query_coords` now documents as its third return value.  So
`cat_from_coords` needed `_, _, IDs = self.query_coords(...)`.  I did not have
to guess.  `cat_from_coords` has been dead for nine years; it now has a test.

**5. SSA.**  The supported route in astropy 8 is `VOTableFile._get_version_checks()`,
which derives all five `version_1_N_or_later` flags from the file's own
`version` (astropy uses it internally at `tree.py:4281`).  It is private, but
hand-rolling the dict is what broke, and adding a second magic key would break
again at astropy 9.  SSA warnings went from 21 to 4.

The `W48 content` finding was right and the fix is satisfying: `content` **is**
the INFO element's text body, but astropy exposes it as an attribute you assign,
not a constructor argument -- passing it to `Info(...)` lands it in `**extra`,
warns, and discards it.  So six INFO elements had been silently losing their
text for every IVOA client.  Added a `_mk_info()` helper that assigns it
properly and stringifies `value` (fixing W08 at the same time).  W03 went away
by giving the five `INPUT:*` params explicit IDs, W47 by giving `TargetName` an
`arraysize`.  The remaining 4 warnings are the VOUnit `Angstrom` deprecation --
see Q&A #18.

**7. Packaging.**  Unlike the `linetools` session, **`build`, `setuptools` and
`setuptools_scm` are all present in `astro`, so I verified by building.**  That
mattered: my first wheel silently contained `specdb/tests/files/tst.xml` and
`specdb/build/tests/files/tst_db.hdf5` -- stale test artifacts picked up by the
package-data glob.  A green test run would never have shown that.  The second
build was clean only after I also cleared setuptools' `build/` scratch
directory, which caches the previous file list.

`setuptools_scm` derives `1.1.dev84+gdef79be47` from the `v1.0` tag, so Q&A #2
is partly answered by the tooling.  `specdb/__init__.py` (previously empty) now
exports `__version__` via `importlib.metadata`, and `docs/conf.py` reads that
instead of its own hard-coded `0.1`.  `*/version.py` added to `.gitignore`.

I deleted `bin/` rather than leave nine shims that are no longer installed --
see Q&A #16.  Package data now includes the test fixtures that `setup.py`'s
`os.walk` omitted, so an installed `specdb` can run its own suite; whether that
is wanted is still Q&A #6.

**8. CI.**  New `.github/workflows/ci_tests.yml`: Python 3.11-3.14 on Linux,
3.13 on macOS, plus a `sphinx-build -W` docs job.  `.travis.yml` deleted.  I
validated the YAML with `yaml.safe_load` and confirmed both branch names parse
as **strings** -- the `'2026'`-as-integer trap from the `linetools` session is
avoided.  CI installs `linetools` from git master, as Travis did, because there
is no compatible release (Q&A #4, #17).

**9. Docs.**  `sphinx-build` now succeeds under `-W`, from a starting point of
not building at all.  Beyond the doc's list I had to fix: a duplicate
`meta-desc` label in `groups.rst`/`meta.rst`; two cross-document `` `name`_ ``
references in `sdss.rst` that docutils cannot resolve (they needed `:ref:`); a
`::` in `casbah.rst` promising a literal block that never arrives (the page is
truncated -- I marked it incomplete); an `ssa.rst` toctree pointing at an
`SSA_IVOA` document that was never generated; and `coords.rst`, which was in no
toctree at all.  Added `.readthedocs.yaml` -- without it RTD cannot rebuild,
which is the actual reason the published docs are frozen at 2022.

**One correction to my own earlier report**, recorded as Q&A #20: the
`download-scripts` label *does* exist, at `docs/scripts.rst:22`.  Only the
reference was misspelled.

**10. Test hygiene.**  `tmp.hdf5`, `tst_db.hdf5` and `tst.xml` all moved to
`tmp_path`.  Verified the way the doc asked: `git status --porcelain` captured
before and after a full run, and the diff is empty.

All three dead assertions were repaired and **all three pass**, so the code
under them was correct and the assertions were simply worthless.  In particular
`assert IDs[2] == IDs[3]` (rows 2 and 3 are duplicate coordinates) passes, so
`get_new_ids` was right all along.

`test_utils.py`'s gate now reads `SPECDB` instead of `IGMSPEC_DB`; I kept the
branch rather than deleting it, because with the real database present it is a
genuine integration test -- it was only ever dead because it watched a variable
nothing sets.  Removed the string-literal-disabled block referencing
`coords_to_spectra`, which exists nowhere.

`.gitignore`: `./build/` became the anchored `/build/`.  I confirmed with
`git check-ignore -v` that `specdb/build/utils.py` is **not** ignored -- a bare
`build/` would have hidden the subpackage.

**11. Python-2 scaffolding.**  25 `from __future__` imports, 10
`# TEST_UNICODE_LITERALS` markers, 3 `basestring` shims, 4 unused `ustr` shims,
2 `bstr` shims, 6 `str('literal')` wrappings, one commented `reload`.  Also
`verbose & (...)` -> `and` in `specdb.py`, and the integer type checks widened
to `np.integer`.  Confirmed the doc's claim that there are no `~`-on-bool bugs
here: all 20 `~` sites are on numpy bool arrays.

**12. Bytes and str.**  Both unambiguous bugs fixed.  The `zsource` truncation
is the one that matters: `np.array(['NONENONE']*n)` is `<U8`, and `BOSS_PIPE`
(a real `z_priority` entry, 9 characters) was being written as `BOSS_PIP` into
every private database built with `mk_meta`.  I first used `dtype=object`, which
**broke `test_mkdb`** -- h5py cannot write object arrays -- so it is `<U32`.
Worth noting as a small lesson: the obvious Python fix was the wrong one for a
package whose output goes to HDF5.

`hdf_decode` no longer raises on a bytes-keyed dict, and now removes the
undecoded key rather than leaving both.

**What I did not do.**  Prompt #6 in full (Q&A #1), and prompt #12 item 1, the
catalog `str`/`bytes` decision (Q&A #7).  Both are behaviour or policy calls
that the document itself says to ask about first.

**No git commands were run** beyond read-only `status`, `log`, `show`,
`check-ignore` and `ls-tree`.  I did delete untracked test artifacts
(`tmp.hdf5`, `tst_db.hdf5`, `tst.xml`) and setuptools' `build/` scratch
directory in order to verify the artifact and packaging work; those were
scratch files, not repository content.
