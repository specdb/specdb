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
        # Name
        self.name = 'SSAI_'+specdb.name

    def querydata(self, POS=None, SIZE=None, TIME=None, BAND=None, FORMAT='HDF5',
                  TOP=None, MAXREC=5000, TARGETCLASS='QSO'):
        """ Perform an SSA-like query on the specdb catalog

        Parameters
        ----------
        POS : str, optional
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

        # METADATA??
        if FORMAT == 'METADATA':
            votable = build_metaquery(def_infos)
            return votable

        def_infos.append(Info(name='REQUEST', value='queryData'))
        def_infos.append(Info(name='serviceName', value='ssap'))
        def_infos.append(Info(name='FORMAT', value=FORMAT))

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

        status = 'OK'
        query_status = []
        qinfos = []

        # Parse POS
        if POS is None:
            status = 'ERROR'
            qinfos.append(Info(name='QUERY_STATUS', value="ERROR", content="POS not provided and not a METADATA query"))
            spl = []
        else:
            spl = POS.split(';')
            scoord = spl[0]
        def_infos.append(Info(name='POS', value=POS))

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


        # Perform query
        _, subcat, IDs = self.specdb.qcat.query_position(coord, SIZE*size_unit, max_match=MAXREC)

        if IDs.size > 0:
            # Grab meta params
            metaparams, pIDs = metaquery_param()

            # Grab sub-catalog
            gd_groups = self.specdb.qcat.groups_containing_IDs(IDs)

            # Loop on groups
            all_vometa = []
            for group in gd_groups:
                # Grab meta from group
                flag_group = self.specdb.qcat.group_dict[group]
                gdID = np.where(subcat['flag_group'].data & flag_group)[0]
                meta_group = self.specdb[group].meta_from_ids(IDs[gdID], first=False)
                meta_attr = self.specdb[group].meta_attr
                # Convert to SSA VO
                ssavo_meta = meta_to_ssa_vo(group, meta_group, meta_attr, subcat[gdID],
                                            self.specdb.idkey,
                                            self.specdb.qcat.cat_attr)
                all_vometa.append(ssavo_meta)
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
                    if hasattr(param, 'unit'):
                        field.unit = param.unit
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
        txt = '<{:s}:  SSA Interface to specdb>'.format(self.name)
        return (txt)

def ssa_defs():
    """ Default dict of SSA compliant spectral info
    Returns
    -------

    """
    ssa_dict = {}
    ssa_dict['DataModel'] = str('Spectrum-2.0')
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


def meta_to_ssa_vo(group, meta, meta_attr, subcat, idkey, cat_attr):
    """ Use a specdb meta table to generate a new astropy Table
    that is ready for conversion to a VOTable.
    Parameters
    ----------
    group : str
      Name of group
    meta : Table
    ssa_attr : dict
    idkey : str

    Returns
    -------
    votbl : astropy Table
      Ready for conversion to VO

    """
    from astropy.time import Time

    # Allow for multiple entries in meta_attr (one per instrument)
    ssa_list = [key for key in meta_attr.keys() if 'SSA_' in key]
    if len(ssa_list) > 0:
        votbls = []
        for key in ssa_list:
            instr = key.split('_')[-1]
            # Cut on Instr
            gd_i = meta['INSTR'] == instr
            if np.sum(gd_i) > 0:
                tdict = {}
                tdict['SSA'] = meta_attr[key].copy()
                votbls.append(meta_to_ssa_vo(group, meta[gd_i], tdict, subcat, idkey, cat_attr))
        votbl = vstack(votbls)
    else:
        ssa_dict = ssa_defs()
        # Get started
        votbl = Table()
        # Dataset
        votbl['DataModel'] = [ssa_dict['DataModel']]*len(meta)
        votbl['DatasetType'] = ssa_dict['DatasetType']
        # DataID
        votbl['Title'] = str(meta_attr['SSA']['Title'])
        votbl['Instrument'] = meta['INSTR']
        # Curation
        votbl['Publisher'] = str(cat_attr['Publisher'])
        # Coord sys
        votbl['SpaceFrameName'] = str(cat_attr['SpaceFrame'])
        votbl['SpaceFrameEquinox'] = cat_attr['EQUINOX']
        # Target
        tnames = []
        for row in meta:
            tnames.append(str('{:s}_{:d}'.format(group, row['GROUP_ID'])))
        votbl['TargetName'] = tnames
        # FluxAxis
        votbl['FluxAxisUcd'] = str(meta_attr['SSA']['FluxUcd'])
        votbl['FluxAxisUnit'] = str(meta_attr['SSA']['FluxUnit'])
        votbl['FluxAxisCalibration'] = str(meta_attr['SSA']['FluxCalib'])
        # SpectraAxis
        votbl['SpectralAxisUcd'] = str(meta_attr['SSA']['SpecUcd'])
        votbl['SpectralAxisUnit'] = str(meta_attr['SSA']['SpecUnit'])
        # Date (MJD)
        dates = Time(meta['DATE-OBS'])
        dates.format = 'mjd'
        votbl['TimeLocation'] = dates.value
        # Pull RA,DEC from main catalog
        radec = np.zeros((len(meta),2))
        for kk,row in enumerate(meta):
            mt = np.where(subcat[idkey] == row[idkey])[0]
            radec[kk,0] = subcat['RA'][mt[0]]
            radec[kk,1] = subcat['DEC'][mt[0]]
        votbl['SpatialLocation'] = radec
        # Spectral location and bounds
        votbl['SpectralLocation'] = (meta['WV_MIN'] + meta['WV_MAX'])/2.
        votbl['SpectralBoundsExtent'] = meta['WV_MAX'] - meta['WV_MIN']
        votbl['SpectralBoundsStart'] = meta['WV_MIN']
        votbl['SpectralBoundsStop'] = meta['WV_MAX']

        # Check against parameters -- Order too
        all_params, pIDs = metaquery_param()
        vo_keys = votbl.keys()
        for vokey,pID in zip(vo_keys,pIDs):
            try:
                assert vokey == pID
            except AssertionError:
                print("{:s} does not match pID".format(vokey))
    # Return
    return votbl


def build_metaquery(def_infos):
    """ Generate a VOTable for a METADATA query
    Parameters
    ----------
    def_infos

    Returns
    -------
    votable : VOTable

    """
    # Begin
    votable = empty_vo(rtype='results')
    resource = votable.resources[0]
    # Info
    sinfo = Info(name='QUERY_STATUS', value="OK", content="Successful search")
    resource.infos.append(sinfo)
    for dinfo in def_infos:
        resource.infos.append(dinfo)
    # Input Params
    iparams = input_params(votable)
    for iparam in iparams:
        resource.params.append(iparam)
    # Output Param
    output_params, _ = metaquery_param(votable)
    for oparam in output_params:
        resource.params.append(oparam)
    # Return
    return votable


def input_params(votbl=None):
    """
    Parameters
    ----------
    votbl

    Returns
    -------

    """
    if votbl is None:
        votbl = empty_vo(rtype='meta')
    all_params = []
    # POS
    pos = Param(votbl, name="INPUT:POS", value="", datatype="char", arraysize="*")
    pos.description = ('The center of the region of interest. The coordinate values are '+
                       'specified in list format (comma separated) in decimal degrees with '+
                       'no embedded white space followed by an optional coord. system.  '+
                       'Allowed systems are (ICRS) and the default is ICRS.')
    all_params.append(pos)
    # SIZE
    size = Param(votbl, name="INPUT:SIZE", value="0.1", datatype="double", unit="deg")
    size.description = ('The radius of the circular region of interest in decimal degrees.'+
                        'Default sized radius is 0.001 degrees')
    all_params.append(size)
    # BAND
    band = Param(votbl, name="INPUT:BAND", value="ALL", datatype="char", arraysize="*")
    band.description = 'Not currently implemented'
    all_params.append(band)
    # TIME
    ptime = Param(votbl, name="INPUT:TIME", value="", datatype="char", arraysize="*")
    ptime.description = 'Not currently implemented'
    all_params.append(ptime)
    # FORMAT
    format = Param(votbl, name="INPUT:FORMAT", value="ALL", datatype="char", arraysize="*")
    format.description = ('Desired format of retrieved data. \n'+
                          'Allowed values are HDF5, METADATA')
    all_params.append(format)
    # Return
    return all_params

def metaquery_param(evotbl=None):
    """
    Parameters
    ----------
    evotbl

    Returns
    -------
    all_params : list
    pIDs : list

    """
    if evotbl is None:
        evotbl = empty_vo(rtype='meta')
    all_params = []
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
    instr = Param(evotbl, ID="Instrument", name="OUTPUT:Instrument",
                  datatype="char", ucd="meta.id;instr", utype="ssa:DataID.Title", arraysize="*", value="")
    instr.description = 'Instrument name'
    all_params.append(instr)

    # data model metadata: Curation.*
    pub_param = Param(evotbl, ID="Publisher", name="OUTPUT:Publisher", utype="ssa:Curation.Publisher", ucd=" meta.curation",
                      datatype="char", arraysize="*", value="")
    all_params.append(pub_param)

    # data model metadata: CoordSys.*
    space_frame = Param(evotbl, ID="SpaceFrameName", name="OUTPUT:SpaceFrameName",
                        datatype="char", utype="ssa:CoordSys.SpaceFrame.Name", arraysize="*", value="")
    space_frame.description = 'Spatial coordinate frame name'
    all_params.append(space_frame)
    equinox = Param(evotbl, ID="SpaceFrameEquinox", name="OUTPUT:SpaceFrameEquinox", datatype="double",
                    ucd="time.equinox;pos.frame", utype="ssa:CoordSys.SpaceFrame.Equinox", unit="yr", value=0.)
    all_params.append(equinox)

    # data model metadata: Target.*
    tname = Param(evotbl, ID="TargetName", name="OUTPUT:TargetName", datatype="char",
                    ucd="meta.id;src", utype="ssa:Target.Name", value="")
    all_params.append(tname)

    # characterization metadata: FluxAxis
    flux_axis_ucd = Param(evotbl, ID="FluxAxisUcd", name="OUTPUT:FluxAxisUcd",
                      datatype="char", utype="ssa:Char.FluxAxis.ucd", arraysize="*", value="")
    all_params.append(flux_axis_ucd)
    flux_axis_unit = Param(evotbl, ID="FluxAxisUnit", name="OUTPUT:FluxAxisUnit",
                          datatype="char", utype="ssa:Char.FluxAxis.unit", arraysize="*", value="")
    all_params.append(flux_axis_unit)
    flux_axis_calib = Param(evotbl, ID="FluxAxisCalibration", name="OUTPUT:FluxAxisCalibration",
                           datatype="char", utype="ssa:Char.FluxAxis.Calibration", arraysize="*", value="")
    all_params.append(flux_axis_calib)

    # characterization metadata: SpectralAxis
    spec_axis_ucd = Param(evotbl, ID="SpectralAxisUcd", name="OUTPUT:SpectralAxisUcd",
                          datatype="char", utype="ssa:Char.SpectralAxis.ucd", arraysize="*", value="")
    all_params.append(spec_axis_ucd)
    spec_axis_unit = Param(evotbl, ID="SpectralAxisUnit", name="OUTPUT:SpectralAxisUnit",
                          datatype="char", utype="ssa:Char.SpectralAxis.unit", arraysize="*", value="")
    all_params.append(spec_axis_unit)

    # characterization metadata: Char.*.Coverage
    time_loc = Param(evotbl, ID="TimeLocation", name="OUTPUT:TimeLocation", datatype="double",
                    ucd="time.epoch", utype="ssa:Char.TimeAxis.Coverage.Location.Value", unit="d", value=0.)
    time_loc.description = "Estimated UT day of observations"
    all_params.append(time_loc)

    sp_loc = Param(evotbl, ID="SpatialLocation", name="OUTPUT:SpatialLocation",
                   datatype="double", ucd="pos.eq", utype="ssa:Char.SpatialAxis.Coverage.Location.Value",
                   arraysize="2", unit="deg", value="", config=cdict)
    sp_loc.description = 'Spatial Position'
    all_params.append(sp_loc)

    spec_loc = Param(evotbl, ID="SpectralLocation", name="OUTPUT:SpectralLocation",
                     datatype="double", ucd="instr.bandpass", utype="ssa:Char.Spectral.Coverage.Location.Value",
                     unit="Angstrom", value=0.)
    spec_loc.description = 'Midpoint of the spectral coverage in Angstroms'
    all_params.append(spec_loc)
    spec_bounds_ext = Param(evotbl, ID="SpectralBoundsExtent", name="OUTPUT:SpectralBoundsExtent",
                   datatype="double", ucd="instr.bandwidth", utype="ssa:Char.Spectral.Coverage.Bounds.Extent",
                   unit="Angstrom", value=0.)
    spec_bounds_ext.description = 'Width of spectrum in Angstroms'
    all_params.append(spec_bounds_ext)
    spec_bounds_strt = Param(evotbl, ID="SpectralBoundsStart", name="OUTPUT:SpectralBoundsStart",
                            datatype="double", ucd="em.wl;stat.min", utype="ssa:Char.Spectral.Coverage.Bounds.Start",
                            unit="Angstrom", value=0.)
    spec_bounds_strt.description = 'Start of spectral coordinate in Angstroms'
    all_params.append(spec_bounds_strt)
    spec_bounds_end = Param(evotbl, ID="SpectralBoundsStop", name="OUTPUT:SpectralBoundsStop",
                             datatype="double", ucd="em.wl;stat.max", utype="ssa:Char.Spectral.Coverage.Bounds.Stop",
                             unit="Angstrom", value=0.)
    spec_bounds_end.description = 'End of spectral coordinate in Angstroms'
    all_params.append(spec_bounds_end)



    # IDs
    pIDs = []
    for param in all_params:
        pIDs.append(param.ID)
    return all_params, pIDs


def default_fields(title, flux='normalized', fxcalib=None):
    """
    Parameters
    ----------
    flux

    Returns
    -------
    def_ssa_dict : dict

    """

    def_ssa_dict = dict( Title=title,
                    SpecUcd='em.wl',
                    SpecUnit='Angstrom',
                    )
    if flux == 'normalized':
        def_ssa_dict['FluxUcd']='arith.ratio;phot.flux.density'
        def_ssa_dict['FluxUnit']=''
        def_ssa_dict['FluxCalib']='NORMALIZED'
    elif flux == 'flambda':
        def_ssa_dict['FluxUcd']='phot.fluDens;em.wl'
        def_ssa_dict['FluxUnit']='erg s**(-1) angstrom**(-1)'
        def_ssa_dict['FluxCalib']='RELATIVE'
    # Over-ride
    if fxcalib is not None:
        def_ssa_dict['FluxCalib'] = fxcalib
    # Return
    return def_ssa_dict
