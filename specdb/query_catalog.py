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
    surveys : list
      List of surveys included in the catalog
    """

    def __init__(self, db_file=None, maximum_ram=10.):
        """
        Returns
        -------

        """
        # Init
        self.cat = None
        self.surveys = None
        # Load catalog
        self.load_cat(db_file)
        # Setup
        self.setup()

    def load_cat(self, db_file):
        """ Open the DB catalog file
        Parameters
        ----------
        db_file : str
          Name of DB file

        Returns
        -------

        """
        import json
        #
        print("Using {:s} for the catalog file".format(db_file))
        hdf = h5py.File(db_file,'r')
        self.cat = Table(hdf['catalog'].value)
        self.db_file = db_file
        # Survey dict
        try: # BACKWARD COMP FOR A BIT
            self.survey_dict = json.loads(hdf['catalog'].attrs['SURVEY_DICT'])
        except:
            self.survey_dict = None
        hdf.close()

    def in_surveys(self, input_surveys, return_list=True):
        """ Return a list of input surveys that are in the DB

        Parameters
        ----------
        in_surveys : list or str
          List of one or more surveys
          If str, converted to list
        surveys : list
          List of surveys to compare against
        return_list : bool, optional
          Return input survey(s) as a list?

        Returns
        -------
        out_surveys : list
          List of overlapping surveys between input and DB

        """
        # Checks
        if isinstance(input_surveys, basestring):
            isurveys = [input_surveys]
        elif isinstance(input_surveys, list):
            isurveys = input_surveys
        else:
            raise IOError("input_surveys must be str or list")
        #
        fsurveys = []
        for isurvey in isurveys:
            if isurvey in self.surveys:
                fsurveys.append(isurvey)
        # Return
        return fsurveys

    def cutid_on_surveys(self, surveys, IGM_IDs):
        """ Find the subset of IGM_IDs within a survey list

        Parameters
        ----------
        surveys : list
          List of surveys to consider, e.g. ['BOSS-DR12', 'SDSS_DR7']
        IGM_IDs : int array

        Returns
        -------

        msk : bool array
          True indicates in survey

        """
        good = np.in1d(self.cat['IGM_ID'], IGM_IDs)
        cut_cat = self.cat[good]
        # Flags
        fs = cut_cat['flag_survey']
        msk = np.array([False]*len(cut_cat))
        for survey in surveys:
            flag = self.survey_dict[survey]
            # In the survey?
            query = (fs % (flag*2)) >= flag
            if np.sum(query) > 0:
                msk[query] = True
        gdIDs = cut_cat['IGM_ID'][msk]
        # Return
        final = np.in1d(IGM_IDs, gdIDs)
        return final

    def match_coord(self, cat_coords, toler=0.5*u.arcsec, verbose=True):
        """ Match an input set of SkyCoords to the catalog within a given radius

        Parameters
        ----------
        coords : SkyCoord
          Single or array
        toler : Angle or Quantity, optional
          Tolerance for a match
        verbose : bool, optional

        Returns
        -------
        indices : bool array
          True = match

        """
        # Checks
        if not isinstance(toler, (Angle, Quantity)):
            raise IOError("Input radius must be an Angle type, e.g. 10.*u.arcsec")
        # Match
        idx, d2d, d3d = match_coordinates_sky(self.coords, cat_coords, nthneighbor=1)
        good = d2d < toler
        # Return
        if verbose:
            print("Your search yielded {:d} matches".format(np.sum(good)))
        return self.cat['IGM_ID'][good]

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
            ID_fg = self.cat['IGM_ID'][close][gdz][izfg]
            ID_bg = self.cat['IGM_ID'][idx[close]][gdz][izfg]
        else:
            pdb.set_trace()
        # Reload
        return ID_fg, ID_bg

    def radial_search(self, inp, radius, verbose=True, private=False):
        """ Search for sources in a radius around the input coord

        Parameters
        ----------
        inp : str or tuple or SkyCoord
          See linetools.utils.radec_to_coord
        toler
        verbose

        Returns
        -------

        """
        # Convert to SkyCoord
        coord = ltu.radec_to_coord(inp)
        # Separation
        sep = coord.separation(self.coords)
        # Match
        good = sep < radius
        # Return
        if verbose:
            print("Your search yielded {:d} match[es]".format(np.sum(good)))
        if private:
            return self.cat['PRIV_ID'][good]
        else:
            return self.cat['IGM_ID'][good]

    def get_cat(self, IGM_IDs):
        """ Grab catalog rows corresponding to the input IDs

        Parameters
        ----------
        IGM_IDs : int array

        Returns
        -------
        rows : Table
          Rows of the catalog

        """
        good = np.in1d(self.cat['IGM_ID'], IGM_IDs)
        return self.cat[good]

    def show_cat(self, IGM_IDs):
        """  Show the catalog

        Parameters
        ----------
        IGM_IDs : int array

        Returns
        -------

        """
        # IGMspec catalog
        good = np.in1d(self.cat['IGM_ID'], IGM_IDs)

        # Catalog keys
        cat_keys = ['IGM_ID', 'RA', 'DEC', 'zem', 'flag_survey']
        for key in self.cat.keys():
            if key not in cat_keys:
                cat_keys += [key]
        self.cat[cat_keys][good].pprint(max_width=120)
        # Print survey dict
        print("----------")
        print("Survey key:")
        for survey in self.surveys:
            print("    {:s}: {:d}".format(survey, self.survey_dict[survey]))
            #print("    {:s}: {:d}".format(survey, idefs.get_survey_dict()[survey]))

    def setup(self):
        """ Set up a few things, e.g. SkyCoord for the catalog
        Returns
        -------

        """
        from specdb import cat_utils as icu
        # SkyCoord
        self.coords = SkyCoord(ra=self.cat['RA'], dec=self.cat['DEC'], unit='deg')
        # Formatting the Table
        self.cat['RA'].format = '8.4f'
        self.cat['DEC'].format = '8.4f'
        self.cat['zem'].format = '6.3f'
        self.cat['sig_zem'].format = '5.3f'
        # Surveys
        unif = np.unique(self.cat['flag_survey'])
        all_surveys = []
        for ifs in unif:
            all_surveys += icu.flag_to_surveys(ifs, self.survey_dict)
        self.surveys = list(np.unique(all_surveys))

    def __repr__(self):
        txt = '<{:s}:  DB_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
