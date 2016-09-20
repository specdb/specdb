""" Module to interface with hdf5 database for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import psutil
import warnings
import h5py
import pdb

import numpy as np

from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D


class InterfaceDB(object):
    """ A Class for interfacing with the DB

    Parameters
    ----------

    Attributes
    ----------
    memory_used : float
      Used memory in Gb
    memory_warning : float
      Value at which a Warning is raised
    hdf : pointer to DB
    maximum_ram : float, optonal
      Maximum memory allowed for the Python session, in Gb
    survey_bool : bool array
      Used to grab data from a given survey
    """

    def __init__(self, db_file, maximum_ram=10.,verbose=True):
        """
        Parameters
        ----------
        db_file : str, optional
        Returns
        -------

        """
        # Init
        self.memory_used = 0.
        self.memory_warning = 5.  # Gb
        self.memory_max = 10.  # Gb
        self.hdf = None
        self.verbose = verbose
        # Open DB
        self.open_db(db_file)
        # Check memory
        self.update()

    def open_db(self, db_file):
        """ Open the DB file
        Parameters
        ----------
        db_file : str

        Returns
        -------

        """
        import json
        #
        if self.verbose:
            print("Using {:s} for the DB file".format(db_file))
        self.hdf = h5py.File(db_file,'r')
        self.db_file = db_file
        self.survey_IDs = None
        #
        self.survey_dict = json.loads(self.hdf['catalog'].attrs['SURVEY_DICT'])
        surveys = self.survey_dict.keys()
        self.surveys = surveys
        if self.verbose:
            print("Available surveys: {}".format(self.surveys))

    def grab_ids(self, survey, IGM_IDs, meta=None, match_meta=None,
                 private=False):
        """ Grab the rows in a survey matching IGM_IDs

        Parameters
        ----------
        survey : str
          Name of the Survey
        IGM_IDs : int or ndarray
          IGM_ID values
        meta : Table, optional
          Meta data for the survey (usually read from hdf)
        match_meta : dict, optional
          key/value to filter spec with, e.g. {'INSTR': HIRES}
        private : bool, optional
          Private DB?

        Returns
        ------
        match_survey : ndarray
          bool array of meta rows from the survey matching IGM_IDs

        """
        # Check
        if survey not in self.hdf.keys():
            raise IOError("Survey {:s} not in your DB file {:s}".format(survey, self.db_file))
        # Find IGMS_IDs indices in survey
        if meta is None:
            meta = Table(self.hdf[survey]['meta'].value)  # This could be too slow..
            meta.meta = dict(survey=survey)
        if private:
            match_survey = np.in1d(meta['PRIV_ID'], IGM_IDs)
        else:
            match_survey = np.in1d(meta['IGM_ID'], IGM_IDs)
        # Match meta?
        if match_meta is not None:
            for key,value in match_meta.items():
                match_survey = match_survey & (meta[key] == value)
        # Store and return
        self.survey_bool = match_survey
        self.meta = meta[match_survey]
        return match_survey

    def grab_meta(self, survey, IGM_IDs=None, show=True):
        """ Grab meta data for survey
        Parameters
        ----------
        survey : str or list
        IGM_IDs : int or array, optional
          Return full table if None
        show : bool, optional
          Show the Meta table (print) in addition to returning

        Returns
        -------
        meta : Table(s)

        """
        #
        if isinstance(survey, list):
            all_meta = []
            for isurvey in survey:
                all_meta.append(self.grab_meta(isurvey, IGM_IDs))
            return all_meta
        # Grab IDs then cut
        meta = Table(self.hdf[survey]['meta'].value)
        if IGM_IDs is not None:
            _ = self.grab_ids(survey, IGM_IDs)
            meta = meta[self.survey_bool]
        else:
            return meta
        self.meta = meta
        self.meta['RA'].format = '7.3f'
        self.meta['DEC'].format = '7.3f'
        self.meta['zem'].format = '6.3f'
        self.meta['WV_MIN'].format = '6.1f'
        self.meta['WV_MAX'].format = '6.1f'

        # Load and return
        return meta

    def grab_spec(self, survey, IGM_IDs, verbose=None, **kwargs):
        """ Grab spectra using staged IDs

        Parameters
        ----------
        survey : str or list
        IGM_IDs : int or intarr

        Returns
        -------
        spec
        meta

        """
        if verbose is None:
            verbose = self.verbose
        if isinstance(survey, list):
            all_spec = []
            all_meta = []
            for isurvey in survey:
                spec, meta = self.grab_spec(isurvey, IGM_IDs, **kwargs)
                if spec is not None:
                    all_spec.append(spec.copy())
                    all_meta.append(meta.copy())
            return all_spec, all_meta
        # Grab IDs
        if self.stage_data(survey, IGM_IDs, **kwargs):
            if np.sum(self.survey_bool) == 0:
                if verbose:
                    print("No spectra matching in survey {:s}".format(survey))
                return None, None
            else:
                if verbose:
                    print("Loaded spectra")
                data = self.hdf[survey]['spec'][self.survey_bool]
        else:
            print("Staging failed..  Not returning spectra")
            return
        # Generate XSpectrum1D
        if 'co' in data.dtype.names:
            co = data['co']
        else:
            co = None
        spec = XSpectrum1D(data['wave'], data['flux'], sig=data['sig'], co=co, masking='edges')
        # Return
        return spec, self.meta

    def show_meta(self, imeta=None, meta_keys=None):
        """ Nicely format and show the meta table
        Parameters
        ----------
        meta_keys : list, optional
          Keys for display
        imeta : Table, optional
          Meta data for the survey (or a subset of it)
          Is pulled from self.meta if not input
        """
        if imeta is None:
            imeta = self.meta
        if meta_keys is None:
            mkeys = ['IGM_ID', 'RA', 'DEC', 'zem', 'SPEC_FILE']
        else:
            mkeys = meta_keys
        #
        for key in imeta.keys():
            if key not in mkeys:
                mkeys += [key]
        imeta[mkeys].pprint(max_width=120)
        return

    def stage_data(self, survey, IGM_IDs, verbose=None, **kwargs):
        """ Stage the spectra for serving
        Mainly checks the memory

        Parameters
        ----------
        survey : str
          Name of the Survey
        IGM_IDs : int or ndarray
          IGM_ID values

        Returns
        -------
        check : bool
        If True, self.survey_IDs is filled
          Indices in the survey dataset

        """
        if verbose is None:
            verbose = self.verbose
        # Checks
        if survey not in self.hdf.keys():
            if survey == 'BOSS_DR12':
                return True
            else:
                raise IOError("Survey {:s} not in your DB file {:s}".format(survey, self.db_file))
        match_survey = self.grab_ids(survey, IGM_IDs, **kwargs)
        # Meta?
        nhits = np.sum(match_survey)
        # Memory check (approximate; ignores meta data)
        spec_Gb = self.hdf[survey]['spec'][0].nbytes/1e9  # Gb
        new_memory = spec_Gb*nhits
        if new_memory + self.memory_used > self.memory_max:
            warnings.warn("This request would exceed your maximum memory limit of {:g} Gb".format(self.memory_max))
            return False
        else:
            if verbose:
                print("Staged {:d} spectra totalling {:g} Gb".format(nhits, new_memory))
            return True


    def update(self):
        """ Update key attributes

        Returns
        -------

        """
        # Memory
        process = psutil.Process(os.getpid())
        self.memory_used = process.memory_info().rss/1e9  # Gb
        if self.memory_used > self.memory_warning:
            warnings.warn("Your memory usage -- {:g} Gb -- is high".format(self.memory_used))

    def __repr__(self):
        txt = '<{:s}:  DB_file={:s} \n'.format(self.__class__.__name__,
                                            self.db_file)
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
