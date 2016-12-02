""" Module to interface with hdf5 database for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import h5py
import numpy as np
import pdb


from astropy.table import Table
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord, match_coordinates_sky, Angle

from linetools import utils as ltu

from specdb.cat_utils import match_ids

try:
    basestring
except NameError:  # For Python 3
    basestring = str


class QueryCatalog(object):
    """ A Class for querying the IGMspec catalog

    Parameters
    ----------

    Attributes
    ----------
    cat : Table
      Astropy Table holding the IGMspec catalog
    groups : list
      List of groups included in the catalog
    """

    def __init__(self, hdf, maximum_ram=10., verbose=False):
        """
        Returns
        -------

        """
        # Init
        self.verbose = verbose
        # Load catalog
        self.load_cat(hdf)
        # Setup
        self.setup()

    def load_cat(self, hdf, idkey=None):
        """ Open the DB catalog file
        Parameters
        ----------
        db_file : str
          Name of DB file
        idkey : str, optional
          Key for ID indices

        Returns
        -------

        """
        import json
        #
        self.cat = Table(hdf['catalog'].value)
        # Set ID key
        self.idkey = idkey
        if idkey is None:
            for key in self.cat.keys():
                if 'ID' in key:
                    if self.idkey is not None:
                        raise ValueError("Two keys with ID in them.  You must specify idkey directly.")
                    self.idkey = key
        # Survey dict
        self.group_dict = json.loads(hdf['catalog'].attrs['SURVEY_DICT'])
        self.groups = list(self.group_dict.keys())
        if self.verbose:
            print("Available groups: {}".format(self.groups))

    def in_groups(self, input_groups, return_list=True):
        """ Return a list of input groups that are in the DB

        Parameters
        ----------
        in_groups : list or str
          List of one or more groups
          If str, converted to list
        groups : list
          List of groups to compare against
        return_list : bool, optional
          Return input group(s) as a list?

        Returns
        -------
        out_groups : list
          List of overlapping groups between input and DB

        """
        # Checks
        if not isinstance(input_groups, basestring):
            igroups = [input_groups]
        elif isinstance(input_groups, list):
            igroups = input_groups
        else:
            raise IOError("input_groups must be str or list")
        #
        fgroups = []
        for igroup in input_groups:
            if igroup in self.groups:
                fgroups.append(igroup)
        # Return
        return fgroups

    def ids_in_groups(self, groups, IDs=None, in_all=False):
        """ Return a list of IDs of sources located in
        one or more groups.  If IDs is input, the subset
        within the groups is returned.

        Default is to require the source occur in at least
        one of the groups.  Use in_all=True to require the
        source be in each of the groups.

        Parameters
        ----------
        groups : list
          List of groups to consider, e.g. ['BOSS-DR12', 'SDSS_DR7']
        IDs : ndarray, optional
          If not input, use the entire catalog of IDs
        in_all : bool, optional
          Require that the source be within *all* of the input groups
          Default is to require it be within at least one group

        Returns
        -------
        gdIDs : int array
        """
        # Init
        ngroup = len(groups)
        if IDs is None:
            IDs = self.cat[self.idkey].data
        # Flags
        cat_rows = match_ids(IDs, self.cat[self.idkey].data)
        fs = self.cat['flag_survey'][cat_rows].data
        msk = np.zeros_like(fs).astype(int)
        for group in groups:
            flag = self.group_dict[group]
            # In the group?
            query = (fs % (flag*2)) >= flag
            msk[query] += 1
        if in_all:
            gd = msk == ngroup
        else:
            gd = msk >= 1
        gdIDs = IDs[gd]
        # Return
        return gdIDs

    def match_coord(self, coords, group=None, toler=0.5*u.arcsec, verbose=True):
        """ Match an input set of SkyCoords to the catalog within a given radius

        Parameters
        ----------
        coords : SkyCoord
          Single or array
        toler : Angle or Quantity, optional
          Tolerance for a match
        group : str, optional
          Restrict to matches within a specific group
        verbose : bool, optional

        Returns
        -------
        indices : int array
          ID values
          -1 if no match within toler
          -2 if within tol but not within input dataset

        """
        # Checks
        if not isinstance(toler, (Angle, Quantity)):
            raise IOError("Input radius must be an Angle type, e.g. 10.*u.arcsec")
        # Match
        idx, d2d, d3d = match_coordinates_sky(coords, self.coords, nthneighbor=1)
        # Generate
        if len(d2d) == 1:
            IDs = np.array([self.cat[self.idkey][idx]])
        else:
            IDs = self.cat[self.idkey][idx].data
        close = d2d < toler
        # Restrict to dataset?
        if group is not None:
            sflag = self.group_dict[group]
            cat_rows = match_ids(IDs, self.cat[self.idkey].data)
            flags = self.cat['flag_survey'][cat_rows]
            query = (flags % (sflag*2)) >= sflag
            IDs[~query] = -2
        # Deal with out of tolerance (after dataset)
        IDs[~close] = -1
        # Finish
        if verbose:
            gd = IDs >= 0
            print("Your search yielded {:d} matches from {:d} input coordinates".format(np.sum(gd),
                                                                                        IDs.size))
        return IDs

    def pairs(self, sep, dv):
        """ Generate a pair catalog
        Parameters
        ----------
        sep : Angle or Quantity
        dv : Quantity
          Offset in velocity.  Positive for projected pairs (i.e. dz > input value)

        Returns
        -------

        """
        # Checks
        if not isinstance(sep, (Angle, Quantity)):
            raise IOError("Input radius must be an Angle type, e.g. 10.*u.arcsec")
        if not isinstance(dv, (Quantity)):
            raise IOError("Input velocity must be a quantity, e.g. u.km/u.s")
        # Match
        idx, d2d, d3d = match_coordinates_sky(self.coords, self.coords, nthneighbor=2)
        close = d2d < sep
        # Cut on redshift
        if dv > 0.:  # Desire projected pairs
            zem1 = self.cat['zem'][close]
            zem2 = self.cat['zem'][idx[close]]
            dv12 = ltu.v_from_z(zem1,zem2)
            gdz = np.abs(dv12) > dv
            # f/g and b/g
            izfg = dv12[gdz] < 0*u.km/u.s
            ID_fg = self.cat[self.idkey][close][gdz][izfg]
            ID_bg = self.cat[self.idkey][idx[close]][gdz][izfg]
        else:
            pdb.set_trace()
        # Reload
        return ID_fg, ID_bg

    def radial_search(self, inp, radius, max=None, verbose=True, private=False):
        """ Search for sources in a radius around the input coord

        Parameters
        ----------
        inp : str or tuple or SkyCoord
          See linetools.utils.radec_to_coord
          Single coordinate
        radius : Angle or Quantity, optional
          Tolerance for a match
        max : int, optional
          Maximum number of matches to return
        verbose

        Returns
        -------
        idx : int array
          Catalog IDs corresponding to match in order of increasing separation
          Returns an empty array if there is no match
        """
        if not isinstance(radius, (Angle, Quantity)):
            raise IOError("Input radius must be an Angle type, e.g. 10.*u.arcsec")
        # Convert to SkyCoord
        coord = ltu.radec_to_coord(inp)
        # Separation
        sep = coord.separation(self.coords)
        # Match
        good = sep < radius
        if verbose:
            print("Your search yielded {:d} match[es] within radius={:g}".format(np.sum(good), radius))
        # Sort by separation
        asort = np.argsort(sep[good])
        if max is not None:
            asort = asort[:max]
        # Return
        return self.cat[self.idkey][good][asort]

    def get_cat(self, IDs):
        """ Grab catalog rows corresponding to the input IDs

        Parameters
        ----------
        IDs : int array

        Returns
        -------
        rows : Table
          Rows of the catalog

        """
        good = np.in1d(self.cat[self.idkey], IDs)
        return self.cat[good]

    def show_cat(self, IDs):
        """  Show the catalog

        Parameters
        ----------
        IDs : int array

        Returns
        -------

        """
        # IGMspec catalog
        good = np.in1d(self.cat[self.idkey], IDs)

        # Catalog keys
        cat_keys = [self.idkey, 'RA', 'DEC', 'zem', 'flag_survey']
        for key in self.cat.keys():
            if key not in cat_keys:
                cat_keys += [key]
        self.cat[cat_keys][good].pprint(max_width=120)
        # Print group dict
        print("----------")
        print("Survey key:")
        for group in self.groups:
            print("    {:s}: {:d}".format(group, self.group_dict[group]))
            #print("    {:s}: {:d}".format(survey, idefs.get_survey_dict()[survey]))

    def setup(self):
        """ Set up a few things, e.g. SkyCoord for the catalog
        Returns
        -------

        """
        # SkyCoord
        self.coords = SkyCoord(ra=self.cat['RA'], dec=self.cat['DEC'], unit='deg')
        # Formatting the Table
        self.cat['RA'].format = '8.4f'
        self.cat['DEC'].format = '8.4f'
        self.cat['zem'].format = '6.3f'
        self.cat['sig_zem'].format = '5.3f'

    def groups_with_IDs(self, IDs, igroup=None):
        """
        Parameters
        ----------
        IDs: int or ndarray
        igroup : list, optional
          List of groups to consider
          Default is the full list of groups

        Returns
        -------
        gd_groups : list
          List of groups containing all of the input IDs

        """
        if isinstance(IDs,int):
            IDs = np.array([IDs])
        nIDs = IDs.size
        if igroup is None:
            igroup = self.groups
        #
        cat_rows = match_ids(IDs, self.cat[self.idkey])
        flags = self.cat['flag_survey'][cat_rows]
        gd_groups = []
        for group in igroup:
            sflag = self.group_dict[group]
            # In the group?
            query = (flags % (sflag*2)) >= sflag
            if np.sum(query) == nIDs:
                gd_groups.append(group)
        # Return
        return gd_groups

    def __repr__(self):
        txt = '<{:s}:  Catalog has {:d} sources\n'.format(self.__class__.__name__,
                                            len(self.cat))
        # Surveys
        txt += '   Loaded groups are {} \n'.format(self.groups)
        txt += '>'
        return (txt)
