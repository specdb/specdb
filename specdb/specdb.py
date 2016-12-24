""" Module for SpecDB Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np
import warnings
import h5py

from astropy import units as u
from astropy.table import Table

from specdb.query_catalog import QueryCatalog
from specdb.interface_group import InterfaceGroup

try:
    basestring
except NameError:  # For Python 3
    basestring = str

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

    def __init__(self, skip_test=True, db_file=None, verbose=False, **kwargs):
        """
        """
        if db_file is None:
            try:
                db_file = self.grab_dbfile(**kwargs)
            except:
                raise IOError("DB not found. Please either check the corresponding environmental "
                              "variable or directly provide the db_file")
        # Init
        self.verbose = verbose
        self.open_db(db_file)
        # Catalog
        self.qcat = QueryCatalog(self.hdf, **kwargs)
        self.cat = self.qcat.cat # For convenience
        self.qcat.verbose = verbose
        self.groups = self.qcat.groups
        self.group_dict = self.qcat.group_dict
        self.idkey = self.qcat.idkey
        # Groups
        self._gdict = {}
        #self.idb = InterfaceDB(db_file, **kwargs)
        #self.idb.idkey = self.qcat.idkey
        #self.idb.verbose = verbose
        # Name, Creation date
        try:
            print("Database is {:s}".format(self.cat.attrs['NAME']))
        except:
            pass
        else:
            print("Created on {:s}".format(self.cat.attrs['CREATION_DATE']))
        # Checks
        '''
        assert self.idb.db_file == self.qcat.db_file
        if not skip_test:
            for group in self.idb.groups:
                try:
                    assert group in self.qcat.groups
                except AssertionError:
                    print("Missing {:s}".format(group))
                    raise IOError
        '''

    def open_db(self, db_file):
        """ Open the DB file

        Parameters
        ----------
        db_file : str

        Returns
        -------

        """
        #
        if self.verbose:
            print("Using {:s} for the DB file".format(db_file))
        self.hdf = h5py.File(db_file,'r')
        self.db_file = db_file

    def ids_to_spectra(self, IDs, group, all_spec=False, chk_group=True, **kwargs):
        """ Grab spectra for input IDs from input dataset
        Default mode is to grab the first spectrum that appears in the dataset
        of a given source (i.e. in cases where more than 1 spectrum exists).
        Use all_spec=True to return a list of XSpectrum1D objects and meta tables,
        one per input source, of all the spectra for each source.

        Parameters
        ----------
        IDs : int ndarray
        group : str
        chk_group : bool, optional
          Check that the IDs are in the group
        all_spec : bool, optional
          Return all spectra for the input coords, un-ordered
        """
        # Check
        if chk_group:
            if not self.qcat.chk_in_group(IDs, group):
                raise IOError("On ore more IDs are not in the input group: {:s}".format(group))
        # Rows
        if all_spec:
            rows = self[group].ids_to_allrows(IDs)
        else:
            rows = self[group].ids_to_firstrow(IDs)
        # Grab and return
        return self[group].grab_specmeta(rows)

    def coords_to_spectra(self, coords, group, tol=0.5*u.arcsec, all_spec=False, **kwargs):
        """ Grab spectra for input coords from input dataset
        Wrapper to ids_to_spectra

        Parameters
        ----------
        coords : SkyCoord
          Typically more than 1
        group : str
        all_spec : bool, optional
          Return all spectra for the input coords, un-ordered
        tol
        kwargs

        Returns
        -------
        spec : XSpectrum1D
          One object containing all of the spectra, in order of input coords
        meta : Table
          Meta data related to spectra
        """
        # Match to catalog
        ids = self.qcat.match_coord(coords, toler=tol, group=group, **kwargs)
        # Check for bad matches
        bad_tol = ids == -1
        if np.sum(bad_tol) > 0:
            print("These input coords have no match in the main catalog")
            print(coords[bad_tol])
            raise IOError("Increase the tolerance for the search or reconsider your query")
        bad_query = ids == -2
        if np.sum(bad_query) > 0:
            print("These input coords are not in the input group {:s}".format(group))
            print(coords[bad_query])
            raise IOError("Try again")
        #
        return self.ids_to_spectra(ids, group, chk_group=False,
                                   all_spec=all_spec, **kwargs)

    def allspec_of_ID(self, ID, groups=None, **kwargs):
        """ Grab all spectra in the database for a given source ID
        Parameters
        ----------
        ID : int
        groups : list, optional
          One or more groups to restrict DB search
        kwargs
          fed to grab_specmeta

        Returns
        -------

        """
        # Restrict groups searched according to user input
        if groups is None:
            groups = self.groups

        # Overlapping groups
        gd_groups = self.qcat.groups_containing_IDs(ID, igroup=groups)

        # Load spectra and meta
        speclist, metalist = [], []
        for group in gd_groups:
            rows = self[group].ids_to_allrows(ID)
            spec, meta = self[group].grab_specmeta(rows, **kwargs)
            # Fill
            speclist.append(spec)
            meta.group = group
            metalist.append(meta)

        # Return
        return speclist, metalist

    def allspec_at_coord(self, coord, tol=0.5*u.arcsec, groups=None, **kwargs):
        """ Radial search for spectra from all data sets for a given coordinate
        Best for single searches (i.e. slower than other approaches)

        Parameters
        ----------
        coord : str, tuple, SkyCoord
          See linetools.utils.radec_to_coord
          Only one coord may be input
        tol : Angle or Quantity, optional
          Search radius
        groups : list, optional
          One or more groups to restrict to
        kwargs :
          fed to grab_specmeta

        Returns
        -------
        speclist : list of XSpectrum1D
          One spectrum per group containing the source
        metalist : list of Tables
          Meta data related to spec

        """
        # Grab ID (will find closest within tol)
        ID = self.qcat.coord_to_ID(coord, tol=tol, **kwargs)

        return self.allspec_of_ID(ID, groups=groups, **kwargs)

    def __getitem__(self, key):
        """ Access the DB groups

        Parameters
        ----------
        key : str

        Returns
        -------

        """
        # Check
        if not isinstance(key, basestring):
            raise IOError("Item must be str")
        # Try to access the dict
        try:
            return self._gdict[key]
        except KeyError:
            if key not in self.groups:
                raise IOError("Input group={:s} is not in the database".format(key))
            else: # Load
                self._gdict[key] = InterfaceGroup(self.hdf, key, idkey=self.idkey)
                return self._gdict[key]

    def __repr__(self):
        txt = '<{:s}:  specDB_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Groups are {} \n'.format(self.groups)
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

    def grab_dbfile(self, version=None, **kwargs):
        """ Grabs the DB file
        Parameters
        ----------
        version : str, optional
          Restrict search to input version

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
        if version is not None:
            fils = glob.glob(db_dir+'/IGMspec_DB_*{:s}*hdf5'.format(version))
        else:
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
        txt += '   Loaded groups are {} \n'.format(self.groups)
        txt += '>'
        return (txt)
