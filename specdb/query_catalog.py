""" Module to interface with hdf5 database for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import h5py
import numpy as np
import pdb
import warnings


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

    def __init__(self, hdf, maximum_ram=10., verbose=False, **kwargs):
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
        # Group dict
        try:
            self.group_dict = json.loads(hdf['catalog'].attrs['GROUP_DICT'])
        except KeyError:
            self.group_dict = json.loads(hdf['catalog'].attrs['SURVEY_DICT'])  # Backwards compatible, will remove
            self.cat['flag_group'] = self.cat['flag_survey']

        self.groups = list(self.group_dict.keys())
        if self.verbose:
            print("Available groups: {}".format(self.groups))

    def chk_in_group(self, IDs, group):
        """ Check whether a set of IDs are in a specified group
        Parameters
        ----------
        IDs : int ndarray
        group : str

        Returns
        -------
        answer : bool
          True if all in group
        in_out : bool ndarray
          True/False for each ID

        """
        # Find rows in catalog
        cat_rows = match_ids(IDs, self.cat[self.idkey].data)
        # Flags
        sflag = self.group_dict[group]
        flags = self.cat['flag_group'][cat_rows]
        # Query on binary
        query = (flags % (sflag*2)) >= sflag
        # Answer
        answer = np.sum(query) == IDs.size
        # Return
        return answer, query

    def coord_to_ID(self, coord, tol=0.5*u.arcsec, closest=True, **kwargs):
        """ Convert an input coord to an ID if matched within a
        given tolerance.  If multiple sources are identified, return
        the closest unless closest=False

        Parameters
        ----------
        coord : str or tuple or SkyCoord
          See linetools.utils.radec_to_coord
          Single coordinate
        tol : Quantity
          Angle
        closest : bool, optional
          If False, raise an error if multiple sources are within tol

        Returns
        -------
        ID : int
          ID of the closest source to the input coord
          within the given tolerance

        """
        # Catalog
        ids = self.radial_search(coord, tol, **kwargs)
        if len(ids) == 0:
            warnings.warn("No sources found at your coordinate within tol={:g}.  Returning None".format(tol))
            return None, None
        elif len(ids) > 1:
            if closest:
                warnings.warn("Found multiple sources in the catalog. Taking the closest one")
            else:
                raise IOError("Multiple sources within tol={:g}.  Refine".format(tol))
        # Finish
        ID = ids[0]
        return ID

    def find_ids_in_groups(self, groups, IDs=None, in_all=True):
        """ Return a list of IDs of sources located in
        one or more groups.  If IDs is input, the subset that are
        within the input groups is returned.

        Default is to require the source occur in all of the input
        the groups.  Use in_all=False to only require the
        source be in at least one of the groups.

        Parameters
        ----------
        groups : list
          List of groups to consider, e.g. ['BOSS-DR12', 'SDSS_DR7']
        IDs : ndarray, optional
          If not input, use the entire catalog of IDs
        in_all : bool, optional
          Require that the source(s) be within *all* of the input groups
          Default is to require it be within at least one group

        Returns
        -------
        gdIDs : int array
          IDs in the group(s)
        good : bool array
          True/False for ID within group(s)
          Mainly useful if user inputs a set of IDs
        """
        # Init
        ngroup = len(groups)
        if IDs is None:
            IDs = self.cat[self.idkey].data
        # Flags
        cat_rows = match_ids(IDs, self.cat[self.idkey].data)
        fs = self.cat['flag_group'][cat_rows].data
        msk = np.zeros_like(fs).astype(int)
        for group in groups:
            flag = self.group_dict[group]
            # In the group?
            query = (fs % (flag*2)) >= flag
            msk[query] += 1
        if in_all:
            good = msk == ngroup
        else:
            good = msk >= 1
        gdIDs = IDs[good]
        # Return
        return gdIDs, good

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
          -2 if within tol but not within input group

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
        # Restrict to group?
        if group is not None:
            answer, query = self.chk_in_group(IDs, group)
            IDs[~query] = -2
        # Deal with out of tolerance (after group)
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

    def radial_search(self, inp, radius, max=None, verbose=True, private=False, **kwargs):
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
        cat_keys = [self.idkey, 'RA', 'DEC', 'zem', 'flag_group']
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

    def groups_containing_IDs(self, IDs, igroup=None):
        """ Return a list of all groups that contain all of the input IDs

        Parameters
        ----------
        IDs: int or ndarray
        igroup : list, optional
          List of groups to consider
          Default is the full list of DB groups

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
        flags = self.cat['flag_group'][cat_rows]
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
