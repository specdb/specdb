""" Module for SpecDB Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np
import warnings

from astropy import units as u
from astropy.table import Table

from specdb.query_catalog import QueryCatalog
from specdb.interface_db import InterfaceDB


class SpecDB(object):
    """ The primary class of this Repository
    Ideally one-stop-shopping for most consumers

    Parameters
    ----------
    skip_test : bool, optional
      Skip tests?  Highly *not* recommended

    Attributes
    ----------
    qcat : QueryCatalog
    idb : InterfaceDB
    """

    def __init__(self, skip_test=True, db_file=None, verbose=None, **kwargs):
        """
        """
        if db_file is None:
            try:
                db_file = self.grab_dbfile()
            except:
                raise IOError("DB not found. Please either check the corresponding environmental "
                              "variable or directly provide the db_file")
        # Init
        self.qcat = QueryCatalog(db_file, **kwargs)
        self.qcat.verbose = verbose
        self.idb = InterfaceDB(db_file, **kwargs)
        self.idb.idkey = self.qcat.idkey
        self.idb.verbose = verbose
        # Name, Creation date
        try:
            print("Database is {:s}".format(self.idb.hdf['catalog'].attrs['NAME']))
        except:
            pass
        else:
            print("Created on {:s}".format(self.idb.hdf['catalog'].attrs['CREATION_DATE']))
        # Checks
        assert self.idb.db_file == self.qcat.db_file
        if not skip_test:
            for survey in self.idb.surveys:
                try:
                    assert survey in self.qcat.surveys
                except AssertionError:
                    print("Missing {:s}".format(survey))
                    raise IOError

    def coords_to_spectra(self, coords, dataset, tol=0.5*u.arcsec, all_spec=False, **kwargs):
        """ Grab spectra for input coords from input dataset
        Default mode is to grab the first spectrum that appears in the dataset
        of a given source (i.e. in cases where more than 1 spectrum exists).
        Use all_spec=True to return a list of XSpectrum1D objects and meta tables,
        one per input source, of all the spectra for each source.


        Parameters
        ----------
        coords : SkyCoord
          Typically more than 1
        dataset : str
        all_spec : bool, optional
        tol
        kwargs

        Returns
        -------
        if all_spec:
        else:
          spec : XSpectrum1D
            One object containing all of the spectra, in order of input coords
          meta : Table
            Meta data related to spectra
        """
        if all_spec:
            meta = Table(self.hdf[dataset]['meta'].value)  # This could be too slow..
            meta.meta = dict(survey=survey)
        else:
            # Match to catalog
            ids = self.qcat.match_coord(coords, tol, dataset=dataset, **kwargs)
            # Check for bad matches
            bad_tol = ids == -1
            if np.sum(bad_tol) > 0:
                print("These input coords have no match in the main catalog")
                print(coords[bad_tol])
                raise IOError("Increase the tolerance for the search or reconsider your query")
            bad_query = ids == -2
            if np.sum(bad_query) > 0:
                print("These input coords are not in the input survey {:s}".format(dataset))
                print(coords[bad_query])
                raise IOError("Try again")
            # Grab and return
            return self.idb.grab_spec(dataset, ids)

    def allspec_at_coord(self, coord, tol=0.5*u.arcsec, isurvey=None, **kwargs):
        """ Radial search for spectra from all data sets for a given coordinate
        Best for single searches (i.e. slower than other approaches)

        Parameters
        ----------
        coord : str, tuple, SkyCoord
          See linetools.utils.radec_to_coord
          Only one coord may be input
        tol : Angle or Quantity, optional
          Search radius
        isurvey : str or list, optional
          One or more surveys to restrict to
        kwargs :
          fed to grab_spec


        Returns
        -------
        speclist : list of XSpectrum1D
          One spectrum per survey containing the source
        metalist : list of Tables
          Meta data related to spec

        """
        # Catalog
        ids = self.qcat.radial_search(coord, tol, **kwargs)
        if len(ids) == 0:
            warnings.warn("No sources found at your coordinate within tol={:g}.  Returning None".format(tol))
            return None, None
        elif len(ids) > 1:
            warnings.warn("Found multiple sources in the catalog. Taking the closest one")
        idv = ids[0]

        # Restrict surveys searched according to user input
        if isurvey is None:
            surveys = self.qcat.surveys
        else:
            surveys = self.qcat.in_surveys(isurvey)

        # Overlapping surveys
        gd_surveys = self.qcat.surveys_with_IDs(idv, isurvey=surveys)

        # Load spectra
        speclist, metalist = self.idb.grab_spec(gd_surveys, idv, **kwargs)
        return speclist, metalist

    def __getattr__(self, k):
        """ Overload attributes using the underlying classes

        Parameters
        ----------
        k

        Returns
        -------

        """
        # Try DB first
        try:
            return getattr(self.idb, k)
        except AttributeError:
            # Try qcat last
            return getattr(self.qcat, k)

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)


class IgmSpec(SpecDB):
    """ Main class for using IGMSpec DB
    See SpecDB for full docs

    Parameters
    ----------
    skip_test : bool, optional
      Skip tests?  Highly *not* recommended

    Attributes
    ----------
    """

    def __init__(self, db_file=None, skip_test=True, **kwargs):
        """
        """
        # db_file
        SpecDB.__init__(self, db_file=db_file, skip_test=skip_test, **kwargs)

    def grab_dbfile(self):
        """ Grabs the DB file
        Returns
        -------
        db_file : str
          full path to the DB file

        """
        import os, glob
        if os.getenv('IGMSPEC_DB') is None:
            warnings.warn('Environmental variable IGMSPEC_DB not set. Assuming this is a test')
            import igmspec
            db_dir = igmspec.__path__[0]+'/tests/files/'
        else:
            db_dir = os.getenv('IGMSPEC_DB')
        #
        fils = glob.glob(db_dir+'/IGMspec_DB_*hdf5')
        fils.sort()
        db_file = fils[-1]  # Should grab the latest
        # Return
        return db_file

    def grab_myers(self):
        """ Returns the Myers quasar catalog as an astropy Table

        Returns
        -------
        myers : Table

        """
        from astropy.table import Table
        return Table(self.idb.hdf['quasars'].value)

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                                                 self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
