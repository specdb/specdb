""" Module for SpecDB Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import pdb
import numpy as np
import warnings
import h5py

from astropy import units as u
from astropy.table import Table, vstack

from specdb import query_catalog as spdb_qc
from specdb import interface_group as spdb_ig
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
        reload(spdb_qc)
        self.verbose = verbose
        self.open_db(db_file)
        # Catalog
        self.qcat = spdb_qc.QueryCatalog(self.hdf, verbose=self.verbose, **kwargs)
        self.cat = self.qcat.cat  # For convenience
        self.qcat.verbose = verbose
        self.groups = self.qcat.groups
        self.group_dict = self.qcat.group_dict
        self.idkey = self.qcat.idkey
        # Groups
        self._gdict = {}
        # Name, Creation date
        self.name = self.qcat.cat_attr['NAME']
        print("Database is {:s}".format(self.name))
        print("Created on {:s}".format(self.qcat.cat_attr['CREATION_DATE']))
        # Return
        return

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

    def meta_from_coords(self, coords, query_dict=None, groups=None,
                               first=True, **kwargs):
        """ Return meta data for an input set of coordinates
        Parameters
        ----------
        coords : SkyCoord
          Expecting an array of coordinates
        query_dict : dict, optional
        groups : list, optional
          If provided, the meta data of the groups are searched in the list order
        first : bool, optional
          Only provide the first entry found for the source
        kwargs

        Returns
        -------
        matches : bool array
          True if the coordinate + query matches in database
        final_meta : masked Table or list
          If first=True (default), the method returns a masked Table
          with each row aligned to the input coordinates.  Entries
          that do not match are fully masked.
          If first=False, this is a list of meta data Tables one per coordinate.
        """
        # Cut down using source catalog
        matches, matched_cat, IDs = self.qcat.query_coords(coords, query_dict=query_dict,
                                                         groups=groups, **kwargs)
        # Setup
        if query_dict is None:
            query_dict = {}

        # Loop on good ones -- slow
        idx = np.where(IDs >= 0)[0]
        if first:
            all_meta = []
        else:
            all_meta = [None]*matches.size
        for iidx in idx:
            query_dict[self.idkey] = [IDs[iidx]]
            # Groups for me
            sub_groups = []
            for group, bit in self.group_dict.items():
                if np.sum(matched_cat['flag_group'][iidx] & bit) > 0:
                    sub_groups.append(group)
            # If groups was input, order by groups
            if groups is not None:
                new_sub = []
                for group in groups:
                    if group in sub_groups:
                        new_sub.append(group)
                # Replace
                sub_groups = new_sub
            # Query
            meta = self.query_meta(query_dict, groups=sub_groups, **kwargs)
            # First?
            if first:
                if meta is None:
                    matches[iidx] = False
                else:
                    all_meta.append(meta[0:1])
            else:
                all_meta[iidx] = meta
        # Finish and pad misses
        if first:
            if len(all_meta) == 0:
                final_meta = None
            else:
                stack = vstack(all_meta)
                final_meta = Table(np.repeat(np.zeros_like(stack[0]), len(IDs)), masked=True)
                final_meta[np.where(matches)] = stack
                # Mask bad rows but fill in IDs
                for row in np.where(~matches)[0]:
                    final_meta.mask[row] = [True]*len(final_meta.mask[row])
                final_meta[self.idkey][np.where(~matches)] = IDs[~matches]
        else:
            final_meta = all_meta


        # Return
        return matches, final_meta

    def meta_from_position(self, inp, radius, query_dict=None, groups=None, **kwargs):
        """  Retrieve meta data for sources around a position on the sky

        Parameters
        ----------
        inp : coordinate in one of several formats
        radius : Angle or Quantity
        query_dict : dict, optional
        groups : list, optional
        kwargs

        Returns
        -------
        all_meta : Table (likely masked)
          All of the meta data for the sources within the region on the sky
          Column 'GROUP' contains group name
          None if there is no match

        """

        # Cut down using source catalog
        matches, sub_cat, IDs = self.qcat.query_position(inp, radius, query_dict=query_dict,
                                                         groups=groups, **kwargs)
        # Add IDs
        if query_dict is None:
            query_dict = {}
        query_dict[self.idkey] = IDs.tolist()

        # Build up groups (to restrict on those that match)
        sub_groups = []
        for group, bit in self.group_dict.items():
            if np.sum(sub_cat['flag_group'] & bit) > 0:
                sub_groups.append(group)
        # Call (restrict on IDs at the least)
        all_meta = self.query_meta(query_dict, groups=sub_groups, **kwargs)
        # Finish
        return all_meta

    def meta_from_ID(self, ID, query_dict=None, **kwargs):
        """ Grabs all meta data for source with a given ID

        Parameters
        ----------
        ID : int
        kwargs : passed to query_meta

        Returns
        -------
        meta : Table
          Meta data table derived from all IDs input
        """
        #
        if query_dict is None:
            query_dict = {}
        query_dict[self.idkey] = [ID]

        # Call (restrict at least on IDs)
        meta = self.query_meta(query_dict, **kwargs)

        # Return
        return meta

    def query_meta(self, qdict, groups=None, **kwargs):
        """ Return all meta data matching the query dict

        Parameters
        ----------
        qdict : dict
          Query dict
        groups : list, optional
        kwargs

        Returns
        -------
        meta : Table
          Stack of all meta data satisfying the query
          Column 'GROUP' contains group name
          None if there is no match
        """
        # Init
        if groups is None:
            groups = self.groups
        # Loop on groups
        all_meta = []
        for group in groups:
            matches, sub_meta, IDs = self[group].query_meta(qdict, **kwargs)
            if len(sub_meta) > 0:
                # Add group
                sub_meta['GROUP'] = str(group)
                # Add RA/DEC?
                # Append -- Not keeping the empty ones
                all_meta.append(sub_meta)
        # Stack
        if len(all_meta) == 0:
            return None
        else:
            return vstack(all_meta)

    def spectra_from_meta(self, meta, debug=False):
        """ Returns one spectrum per row in the input meta data table
        This meta data table should have been generated by a meta query

        Parameters
        ----------
        meta : Table
          Must include a column 'GROUP' indicating the group for each row
          This automatically generated by any of the meta query methods

        Returns
        -------
        spec : XSpectrum1D
          An object containing all of the spectra
          These are aligned with the input meta table
        """
        from linetools.spectra import utils as ltsu
        # Checks
        if 'GROUP' not in meta.keys():
            print("Input meta data table must include a GROUP column")
            print("We suspect yours was not made by a meta query.  Try one")
            raise IOError("And then try again.")
        #
        groups = np.unique(meta['GROUP'].data)
        all_spec = []
        sv_rows = []
        for group in groups:
            sub_meta = meta['GROUP'] == group
            sv_rows.append(np.where(sub_meta)[0])
            # Grab
            all_spec.append(self[group].spec_from_meta(meta[sub_meta]))
        # Collate
        spec = ltsu.collate(all_spec)
        # Re-order
        idx = np.concatenate(sv_rows)
        srt = np.argsort(idx)
        spec2 = spec[srt]
        #
        return spec2

    def spectra_from_coord(self, inp, tol=0.5*u.arcsec, **kwargs):
        """ Return spectra and meta data from an input coordinate
        Radial search at that location within a small tolerance
        Returns closet source if multiple are found

        Parameters
        ----------
        inp : str, tuple, SkyCoord
          See linetools.utils.radec_to_coord
          Only one coord may be input
        tol : Angle or Quantity, optional
          Search radius
        kwargs :
          fed to grab_specmeta

        Returns
        -------
        spec : XSpectrum1D
          All spectra corresponding to the closest source within tolerance
        meta : Table
          Meta data describing to spec

        """
        # Grab meta
        meta = self.meta_from_position(inp, tol, max_match=1, **kwargs)
        if meta is None:
            return None, None
        # Grab spec and return
        return self.spectra_from_meta(meta), meta

    def spectra_from_ID(self, ID, **kwargs):
        """ Return all spectra for a given source ID

        Parameters
        ----------
        ID : int
        kwargs

        Returns
        -------
        spec : XSpectrum1D
          All spectra corresponding to the closest source within tolerance
        meta : Table
          Meta data describing to spec

        """
        # Grab meta
        meta = self.meta_from_ID(ID, **kwargs)
        # Grab spec and return
        return self.spectra_from_meta(meta), meta

    def spectra_in_group(self, coords, group, **kwargs):
        """ Grab the first spectrum found in a given group for an input set of coordinates
        Parameters
        ----------
        coords : SkyCoord
          Expected to be an array
        group : str
          Name of group to use
        kwargs

        Returns
        -------
        spec : XSpectrum1D
          First spectrum found for each input coordinate.
          Aligned to set of coordinates
        meta : Table
          Meta data describing the spectra; also aligned
        """
        # Grab meta
        matches, meta = self.meta_from_coords(coords, groups=[group], **kwargs)
        # Check
        if np.sum(matches) != matches.size:
            raise IOError("Not all of the input coordinates matched in your group")
        # Grab spectra
        spec = self.spectra_from_meta(meta)
        # Return
        return spec, meta

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
                reload(spdb_ig)
                self._gdict[key] = spdb_ig.InterfaceGroup(self.hdf, key, idkey=self.idkey)
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
