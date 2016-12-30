""" Module for IVOA SSA
See http://www.ivoa.net/documents/SSA/
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import h5py
import numpy as np
import pdb
import warnings

from astropy.io.votable.tree import VOTableFile, Resource, Info, Param
from astropy.io.votable import from_table

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack

from linetools import utils as ltu

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

        if IDs.size > 0:
            # Grab meta params
            evotbl = empty_vo(rtype='meta')
            metaparams, pdict, pIDs = self.metaquery_response(evotbl)

            # Grab sub-catalog
            subcat = self.specdb.qcat.cat_from_ids(IDs)
            gd_groups = self.specdb.qcat.groups_containing_IDs(IDs)

            # Loop on groups
            all_vometa = []
            for group in gd_groups:
                # Grab meta from group
                flag_group = self.specdb.qcat.group_dict[group]
                gdID = np.where(subcat['flag_group'].data & flag_group)[0]
                meta_group = self.specdb[group].cut_meta(IDs[gdID], first=False)
                meta_attr = self.specdb[group].meta_attr
                meta_attr['SSA'] = None
                vo_meta = meta_to_ssa_vo(meta_group, meta_attr['SSA'], subcat[gdID])
                all_vometa.append(vo_meta)
            vometa = vstack(all_vometa)
            # Generate true VOTable
            votable = from_table(vometa)
            # Update fields
            tbl = votable.resources[0].tables[0]
            for param in metaparams:
                try:
                    field = tbl.get_field_by_id(param.ID)
                except KeyError:
                    print("Need field with ID={:s}".format(param.ID))
                else:
                    field.utype = param.utype
                    if hasattr(param, 'ucd'):
                        field.ucd = param.ucd
            # Add Parameters
            #pub_param = Param(tbl, name="Publisher", utype="ssa:Curation.Publisher", ucd=" meta.curation",
            #                  datatype="char", arraysize="*", value="JXP")
            #tbl.params.append(pub_param)
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

def ssa_defs():
    """ Default dict of SSA compliant spectral info
    Returns
    -------

    """
    ssa_dict = {}
    ssa_dict['DataModel'] = str('Spectrum-1.0')
    ssa_dict['DatasetType'] = str('Spectrum')

    return ssa_dict


def empty_vo(rtype='results'):
    """ Generate an empty VO Table with one resource

    Parameters
    ----------
    rtype : str, optional
      Type of resouuce

    Returns
    -------
    evotable : VOTable

    """
    evotable = VOTableFile()
    resource = Resource(type=rtype)
    evotable.resources.append(resource)
    return evotable


def meta_to_ssa_vo(meta, ssa_attr, subcat):
    """
    Parameters
    ----------
    meta : Table
    ssa_attr : dict

    Returns
    -------
    votbl : astropy Table
      Ready for conversion to VO

    """
    ssa_dict = ssa_defs()
    # Get started
    votbl = Table()
    votbl['DataModel'] = [ssa_dict['DataModel']]*len(meta)
    votbl['DatasetType'] = ssa_dict['DatasetType']
    #
    # RA,DEC
    radec = np.zeros((len(meta),2))
    radec[:,0] = subcat['RA']
    radec[:,1] = subcat['DEC']
    votbl['SpatialLocation'] = radec

    # Check against parameters -- Order too
    all_params, param_dict, pIDs = metaquery_param()
    vo_keys = votbl.keys()
    for vokey,pID in zip(vo_keys,pIDs):
        assert vokey == pID
    # Return
    return votbl

def metaquery_param(evotbl=None):
    if evotbl is None:
        evotbl = empty_vo(rtype='meta')
    all_params = []
    param_dict = {}  # Converts ID to specdb tag
    cdict = {'version_1_3_or_later': True}
    # Service metadata

    # data model metadata: Dataset.*
    datamodel = Param(evotbl, ID="DataModel", name="OUTPUT:DataModel", datatype="char",
                      utype="ssa:Dataset.DataModel", arraysize="*", value="")
    datamodel.description = 'Datamodel name and version'
    all_params.append(datamodel)
    dataset = Param(evotbl, ID="DatasetType", datatype="char", name="OUTPUT:DatasetType",
                    utype="ssa:Dataset.Type", value="Spectrum", arraysize="*")
    dataset.description = 'Dataset type'
    all_params.append(dataset)

    # data model metadata: DataID.*
    title = Param(evotbl, ID="Title", name="OUTPUT:Title",
                  datatype="char", ucd="meta.title;meta.dataset", utype="ssa:DataID.Title", arraysize="*", value="")
    title.description = 'Dataset Title'
    all_params.append(title)

    # characterization metadata: Char.*.Coverage
    sp_loc = Param(evotbl, ID="SpatialLocation", name="OUTPUT:SpatialLocation", datatype="double", ucd="pos.eq", utype="ssa:Char.SpatialAxis.Coverage.Location.Value", arraysize="2", unit="deg", value="", config=cdict)
    sp_loc.description = 'Spatial Position'
    all_params.append(sp_loc)
    param_dict['SpatialLocation'] = ['RA','DEC']


    pub_param = Param(evotbl, ID="Publisher", name="OUTPUT:Publisher", utype="ssa:Curation.Publisher", ucd=" meta.curation",
                      datatype="char", arraysize="*", value="")
    # IDs
    pIDs = []
    for param in all_params:
        pIDs.append(param.ID)
    return all_params, param_dict, pIDs
