""" Module for IVOA SSA
See http://www.ivoa.net/documents/SSA/
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


class SSAInterface(object):
    """ A Class for handling Simple Spectral Access (SSA) activities

    Parameters
    ----------
    specdb : SpecDB object
    """

    def __init__(self, specdb, maximum_ram=10., verbose=False, **kwargs):
        """
        Returns
        -------

        """
        # Init
        self.verbose = verbose
        # Pointer
        self.specdb = specdb

    def querydata(self, POS, SIZE=None, TIME=None, BAND=None, FORMAT='HDF5',
                  TOP=None, MAXREC=5000, TARGETCLASS='QSO'):
        """ Perform an SSA-like query on the specdb catalog

        Parameters
        ----------
        POS : str
          position on the sky RA,DEC;coordinate system
          Only ICRS is accepted for now
        SIZE : float, optional
          Search radius in deg
        TIME : str, optional
          Not implemented currently
        BAND : str, optional
        FORMAT : str, optional
          Specifies format of dataset that would be returned

        Returns
        -------
        result : VOTable

        """
        from astropy.io.votable.tree import VOTableFile, Resource, Info, Param
        from astropy.io.votable import from_table

        def empty_vo():
            evotable = VOTableFile()
            resource = Resource(type="results")
            evotable.resources.append(resource)
            return evotable

        # Default Infos
        def_infos = []
        def_infos.append(Info(name='SERVICE_PROTOCOL', value=1.1, content="SSAP"))
        def_infos.append(Info(name='REQUEST', value='queryData'))
        def_infos.append(Info(name='serviceName', value='ssap'))
        def_infos.append(Info(name='POS', value=POS))
        def_infos.append(Info(name='FORMAT', value=FORMAT))

        status = 'OK'
        query_status = []
        qinfos = []
        # Parse POS
        spl = POS.split(';')
        scoord = spl[0]

        if len(spl) > 1:
            coord_sys = spl[1]
        else:
            coord_sys = 'ICRS'
        if coord_sys not in ['ICRS']:
            status = 'ERROR'
            qinfos.append(Info(name='QUERY_STATUS', value="ERROR", content="Coordinate system {:s} not implemented".format(coord_sys)))
            #raise IOError("Coordinate system {:s} not implemented".format(coord_sys))

        # Return if failed
        if status != 'OK':
            votable = empty_vo()
            votable.resources[0].infos += qinfos
            votable.resources[0].infos += def_infos
            return votable
        else:
            qinfos.append(Info(name='QUERY_STATUS', value="OK", content="Successful search"))

        # Build Table
        ra,dec = scoord.split(',')
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')

        # SIZE
        size_unit = u.deg
        if SIZE is None:
            SIZE = 1e-3  # deg
        def_infos.append(Info(name='SIZE', value=SIZE, content="Search radius adopted (deg)"))

        # Other
        if TIME is not None:
            warnings.warn("TIME parameter is not yet implemented")
        if BAND is not None:
            warnings.warn("BAND parameter is not yet implemented")

        # Perform query
        IDs = self.specdb.qcat.radial_search(coord, SIZE*size_unit, mt_max=MAXREC)

        #
        if IDs.size > 0:
            # Grab sub-catalog
            subcat = self.specdb.qcat.cat_from_ids(IDs)
            # Add Units
            subcat['RA'].unit = u.deg
            subcat['DEC'].unit = u.deg
            #
            votable = from_table(subcat)
            tbl = votable.resources[0].tables[0]
            # Add description text to Fields
            tbl.get_field_by_id("RA").description = 'Right Ascension (J2000)'
            # Add Parameters
            pub_param = Param(tbl, name="Publisher", utype="ssa:Curation.Publisher", ucd=" meta.curation",
                              datatype="char", arraysize="*", value="JXP")
            tbl.params.append(pub_param)
        else:  # Generate a dummy table
            votable = empty_vo()

        # INFO
        for qinfo in qinfos:
            votable.resources[0].infos.append(qinfo)
        for dinfo in def_infos:
            votable.resources[0].infos.append(dinfo)

        # Return
        return votable


    def __repr__(self):
        txt = '<{:s}:  SSA Interface to specdb>'
        return (txt)
