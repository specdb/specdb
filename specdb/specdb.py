""" Module for SpecDB Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

#import numpy as np
import pdb
import warnings

from astropy import units as u

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

    def spec_from_coord(self, coord, tol=5.*u.arcsec, isurvey=None, **kwargs):
        """ Radial search for spectra around given coordinate
        Best for single searches (i.e. slower than other approaches)

        Parameters
        ----------
        coord : str, tuple, SkyCoord
          See linetools.utils.radec_to_coord
        tol : Angle or Quantity, optional
          Search radius
        isurvey : str or list, optional
          One or more surveys to include
        kwargs :
          fed to grab_spec


        Returns
        -------
        spec : XSpectrum1D or list of XSpectrum1D
          One or more spectra satisfying the radial search
        meta : Table or list of Tables
          Meta data related to spec

        """
        # Catalog
        ids = self.qcat.radial_search(coord, tol, **kwargs)
        if len(ids) == 0:
            warnings.warn("No sources found at your coordinate.  Returning none")
            return None, None
        elif len(ids) > 1:
            warnings.warn("Found multiple sources.  Hope you expected that.")

        # Surveys
        if isurvey is None:
            surveys = self.qcat.surveys
        else:
            surveys = self.qcat.in_surveys(isurvey)

        # Load spectra
        spec, meta = self.idb.grab_spec(surveys, ids, **kwargs)
        return spec, meta

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

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                                                 self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
