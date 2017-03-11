""" Module to interface with an hdf5 data group in specdb
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

from specdb.cat_utils import match_ids
from specdb.group_utils import show_group_meta
from specdb import utils as spdbu

class InterfaceGroup(object):
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
    """

    def __init__(self, hdf, group, idkey, maximum_ram=10., verbose=True, **kwargs):
        """
        Parameters
        ----------
        hdf : h5py.File object
        group : str
        idkey : str

        Returns
        -------

        """
        # Init
        self.hdf = hdf
        self.group = group
        self.idkey = idkey
        self.verbose = verbose
        # Load meta
        self.load_meta(group, **kwargs)
        # Memory
        self.memory_used = 0.
        self.memory_warning = 5.  # Gb
        self.memory_max = 10.  # Gb
        self.update()

    def load_meta(self, group, reformat=True):
        """ Load the meta data as a Table
        Parameters
        ----------
        group : str
        """
        import json
        self.meta = spdbu.hdf_decode(self.hdf[group+'/meta'].value, itype='Table')
        # Attributes
        self.meta_attr = {}
        for key in self.hdf[group+'/meta'].attrs.keys():
            if 'SSA' in key:
                self.meta_attr[key] = json.loads(spdbu.hdf_decode(self.hdf[group+'/meta'].attrs[key]))
            else:
                self.meta_attr[key] = spdbu.hdf_decode(self.hdf[group+'/meta'].attrs[key])
        # Reformat
        if reformat:
            try:
                self.meta['RA_GROUP'].format = '8.4f'
            except KeyError:  # Backwards compatible, will deprecate
                self.meta['RA'].format = '8.4f'
                self.meta['DEC'].format = '8.4f'
                self.meta['zem'].format = '6.3f'
            else:
                self.meta['DEC_GROUP'].format = '8.4f'
                self.meta['zem_GROUP'].format = '6.3f'
            self.meta['WV_MIN'].format = '6.1f'
            self.meta['WV_MAX'].format = '6.1f'
        # Add group
        self.meta.meta['group'] = group

    def groupids_to_rows(self, group_IDs):
        """ Convert GROUP_ID values to rows in the meta table
        Mainly used to then grab the corresponding spectra

        Parameters
        ----------
        group_IDs : int or ndarray

        Returns
        -------
        rows : ndarray

        """
        # Checking
        if isinstance(group_IDs, int):
            group_IDs = np.array([group_IDs])  # Insures meta and other arrays are proper
        # Find rows
        rows = match_ids(group_IDs, self.meta['GROUP_ID'])
        # Return
        return rows

    def ids_to_firstrow(self, IDs):
        """ Given an input set of IDs, pass back the array of rows that
        match.  If a given source has more than one entry, only the first
        one is returned.

        Parameters
        ----------
        IDs : int or ndarray
          ID values
          Converted to array if int

        Returns
        -------
        rows : int ndarray
          Array of indices in meta table, that correspond to the input IDs (in order)
        """
        # Checking
        if isinstance(IDs, int):
            IDs = np.array([IDs])  # Insures meta and other arrays are proper
        # Check that input IDs are all covered
        chk_survey = np.in1d(IDs, self.meta[self.idkey])
        if np.sum(chk_survey) != IDs.size:
            raise IOError("Not all of the input IDs are located in requested group: {:s}".format(self.group))
        # All matching rows
        match_group = np.in1d(self.meta[self.idkey], IDs)
        # Find indices of input IDs in meta table -- first instance in meta only!
        gdi = self.meta[self.idkey].data[match_group]
        xsorted = np.argsort(gdi)
        ypos = np.searchsorted(gdi, IDs, sorter=xsorted)
        indices = xsorted[ypos]  # Location in subset of meta table
        # Store and return
        rows = np.where(match_group)[0][indices]
        return rows

    def ids_to_allrows(self, IDs):
        """ Identify all the rows in a dataset matching input IDs
        All of the rows matching the input IDs are returned
        The current implementation checks that all sources exist
        within the data group.

        Parameters
        ----------
        IDs : int or ndarray
          ID values
          Converted to array if int

        Returns
        ------
        rows : int ndarray
          Array of all rows in meta table that match to the input IDs
          Ordered as in meta table
          This array will exceed the size of the input array if there
          is more than one spectrum per source
        """
        # Checking
        if isinstance(IDs, int):
            IDs = np.array([IDs])  # Insures meta and other arrays are proper
        # Check that input IDs are all covered
        try:
            chk_survey = np.in1d(IDs, self.meta[self.idkey])
        except:
            pdb.set_trace()
        if np.sum(chk_survey) != IDs.size:
            raise IOError("Not all of the input IDs are located in requested group: {:s}".format(self.group))
        # Find rows (bool array)
        match_survey = np.in1d(self.meta[self.idkey], IDs)
        return np.where(match_survey)[0]

    def grab_specmeta(self, rows, verbose=None, **kwargs):
        """ Grab the spectra and meta data for an input set of rows
        Aligned to the rows input

        Parameters
        ----------
        rows : int or ndarray
        verbose
        kwargs

        Returns
        -------
        spec : XSpectrum1D
          Spectra requested, ordered by the input rows
        meta : Table  -- THIS MAY BE DEPRECATED
          Meta table, ordered by the input rows
        """
        if isinstance(rows, int):
            rows = np.array([rows])  # Insures meta and other arrays are proper
        if verbose is None:
            verbose = self.verbose
        # Check memory
        if self.stage_data(rows, **kwargs):
            if verbose:
                print("Loaded spectra")
            # Load
            msk = np.array([False]*len(self.meta))
            msk[rows] = True
            tmp_data = self.hdf[self.group]['spec'][msk]
            # Replicate and sort according to input rows
            idx = match_ids(rows, np.where(msk)[0])
            data = tmp_data[idx]
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
        return spec, self.meta[rows]

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

    def meta_from_ids(self, IDs, first=True):
        """ Grab meta data given a list of IDs

        Parameters
        ----------
        IDs : int or array
          Return full table if None
        first : bool, optional
          Grab only the first row matching?

        Returns
        -------
        meta : Table(s)
          If first=True, aligned to input IDs

        """
        # Grab meta table
        if first:
            rows = self.ids_to_firstrow(IDs)
        else:
            rows = self.ids_to_allrows(IDs)
        cut_meta = self.meta[rows]
        return cut_meta

    def query_meta(self, qdict, **kwargs):
        """
        Parameters
        ----------
        qdict : dict
          Query_dict

        Returns
        -------
        matches : bool array
        sub_meta : Table
          Subset of the meta table matching the query
        IDs : int ndarray
        """
        # Query
        matches = spdbu.query_table(self.meta, qdict, tbl_name='meta data')

        # Return
        return matches, self.meta[matches], self.meta[self.idkey][matches].data

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
        # Show
        show_group_meta()

    def spec_from_meta(self, meta):
        """ Return spectra aligned to input meta table

        Parameters
        ----------
        meta : Table

        Returns
        -------
        spectra : XSpectrum1D

        """
        rows = self.groupids_to_rows(meta['GROUP_ID'])
        # Grab spectra
        spec, _ = self.grab_specmeta(rows)
        # Return
        return spec

    def stage_data(self, rows, verbose=None, **kwargs):
        """ Stage the spectra for serving
        Mainly checks the memory

        Parameters
        ----------
        rows : ndarray
          Indices of desired data

        Returns
        -------
        check : bool
        If True, self.survey_IDs is filled
          Indices in the survey dataset

        """
        if verbose is None:
            verbose = self.verbose
        # Memory check (approximate; ignores meta data)
        spec_Gb = self.hdf[self.group]['spec'][0].nbytes/1e9  # Gb
        new_memory = spec_Gb*rows.size
        if new_memory + self.memory_used > self.memory_max:
            warnings.warn("This request would exceed your maximum memory limit of {:g} Gb".format(self.memory_max))
            return False
        else:
            if verbose:
                print("Staged {:d} spectra totalling {:g} Gb".format(len(rows), new_memory))
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
        txt = '<{:s}:  Group={:s} \n'.format(self.__class__.__name__,
                                            self.group)
        return (txt)
