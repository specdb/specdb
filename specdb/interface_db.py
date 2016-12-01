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

    def meta_rows_from_ids(self, meta, IDs):
        """ Grab the rows in a dataset matching input IDs
        All IDs *must* occur within the survey
        All of the rows matching the input IDs are returned

        Parameters
        ----------
        meta : Table
          Meta data for the dataset (taken previously from hdf DB file)
        IDs : int or ndarray
          ID values
          Converted to array if int

        Returns
        ------
        match_survey : ndarray
          bool array of meta rows from the survey matching IDs
        indices :
        Also fills self.meta and self.indices and self.survey_bool

        """
        # For checking
        if isinstance(IDs, int):
            IDs = np.array([IDs])  # Insures meta and other arrays are proper
        # Check that input IDs are all covered
        chk_survey = np.in1d(IDs, meta[self.idkey])
        if np.sum(chk_survey) != IDs.size:
            raise IOError("Not all of the input IDs are located in your survey")
        # Find rows to grab (bool array)
        match_survey = np.in1d(meta[self.idkey], IDs)
        # Find indices of input IDs in meta table -- first instance in meta only!
        gdi = meta[self.idkey].data[match_survey]
        xsorted = np.argsort(gdi)
        ypos = np.searchsorted(gdi, IDs, sorter=xsorted)
        indices = xsorted[ypos] # Location in subset of meta table (spectra to be grabbed) of the input IDs
        # Store and return
        self.survey_bool = match_survey
        self.meta = meta[match_survey][indices] # only the first entry in table
        self.indices = indices
        return match_survey

    def grab_meta(self, survey, IDs=None, reformat=True):
        """ Grab meta data for survey
        Parameters
        ----------
        survey : str
        IDs : int or array, optional
          Return full table if None

        Returns
        -------
        meta : Table(s)

        """
        # Grab meta table
        meta = Table(self.hdf[survey]['meta'].value)
        # Cut on IDs?
        if IDs is not None:
            _ = self.meta_rows_from_ids(meta, IDs)
            meta = meta[self.survey_bool]
        self.meta = meta
        if reformat:
            self.meta['RA'].format = '7.3f'
            self.meta['DEC'].format = '7.3f'
            self.meta['zem'].format = '6.3f'
            self.meta['WV_MIN'].format = '6.1f'
            self.meta['WV_MAX'].format = '6.1f'

        return meta

    def grab_spec(self, survey, IDs, verbose=None, **kwargs):

    def loop_grab_spec(self, survey, IDs, verbose=None, **kwargs):
        """ Grab spectra using staged IDs
        All IDs must occur in each of the surveys listed

        Order of spectra and meta tables will match the input IDs

        Parameters
        ----------
        survey : str or list
        IDs : int or intarr

        Returns
        -------
        spec : XSpectrum1D
          Spectra requested, ordered by the IDs
        meta : Table
          Meta table, ordered by the IDs

        """
        if verbose is None:
            verbose = self.verbose
        if isinstance(survey, list):
            all_spec = []
            all_meta = []
            for isurvey in survey:
                spec, meta = self.grab_spec(isurvey, IDs, **kwargs)
                if spec is not None:
                    all_spec.append(spec.copy())
                    all_meta.append(meta.copy())
            return all_spec, all_meta
        # Grab IDs
        if self.stage_data(survey, IDs, **kwargs):
            if np.sum(self.survey_bool) == 0:
                if verbose:
                    print("No spectra matching in survey {:s}".format(survey))
                return None, None
            else:
                if verbose:
                    print("Loaded spectra")
                tmp_data = self.hdf[survey]['spec'][self.survey_bool]
                # Replicate and sort according to input IDs
                data = tmp_data[self.indices]
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
            mkeys = [self.idkey, 'RA', 'DEC', 'zem', 'SPEC_FILE']
        else:
            mkeys = meta_keys
        #
        for key in imeta.keys():
            if key not in mkeys:
                mkeys += [key]
        imeta[mkeys].pprint(max_width=120)
        return

    def stage_data(self, survey, IDs, verbose=None, **kwargs):
        """ Stage the spectra for serving
        Mainly checks the memory

        Parameters
        ----------
        survey : str
          Name of the Survey
        IDs : int or ndarray
          ID values

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
        match_survey = self.grab_ids(survey, IDs, **kwargs)
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
