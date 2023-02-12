import xarray as xr
import numpy as np
import sys

from CASutils import calendar_utils as cal

avog = 6.022e23 # Avogadro's number
re = 6.3712e8 # Radius of the earth in cm

def convert_molecules_to_tg(dat,varname,molecular_weight):
    """ Convert surface emissions in molecules/cm2/s to Tg 

    """
    # Start with moleculre/cm2/s.  Convert from molecules to grams.  
    # Divide by Avogadro's number to convert from molecules to moles.  
    #Multiply by molcular weight in g/mol to end up with g/cm2/s
    dat_g = molecular_weight*dat[varname]/avog

    # convert from per s to per year
    dat_g_y = dat_g*365.*86400.

    if "altitude" in dat.dims:
        print('you have altitudes')
        dz = dat.altitude_int[1:dat.altitude_int.size].values - dat.altitude_int[0:dat.altitude_int.size-1].values
        dz = xr.DataArray(dz, coords=[dat.altitude], dims=['altitude'], name='dz')
        dz = dz*1000.*100. # convert to cm
        dat_g_y = (dat_g_y*dz).sum('altitude')


    # Integrate over space
    dlon = np.deg2rad( (dat.lon[2] - dat.lon[1]))
    dlat = np.deg2rad( (dat.lat[2] - dat.lat[1]))
    area = xr.ones_like(dat.isel(time=0))
    weights = np.cos(np.deg2rad(area.lat))*dlat*dlon*re**2. # area in cm2
    dat_g_y_w = dat_g_y.weighted(weights)
    dattot = dat_g_y_w.sum(("lon","lat"))

    # Convert from grams to terra grams
    dattot = dattot/1e12

    return dattot


basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM1/piControl/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/"

so4_a1 = xr.open_dataset(basepath+'ext/ar5_mam3_so4_a1_elev_1850_c090726.nc')
time = cal.MMDD2date(so4_a1.date)
molecular_weight = 115.
ene_tg = convert_molecules_to_tg(so4_a1,'emiss_ene', molecular_weight)
ind_tg = convert_molecules_to_tg(so4_a1,'emiss_ind', molecular_weight)
forestfire_tg = convert_molecules_to_tg(so4_a1,'forestfire', molecular_weight)
grassfire_tg = convert_molecules_to_tg(so4_a1,'grassfire', molecular_weight)
contvolc_a1_tg = convert_molecules_to_tg(so4_a1,'contvolc', molecular_weight)
ene_tg['time'] = time
ind_tg['time'] = time
forestfire_tg['time'] = time
grassfire_tg['time'] = time
contvolc_a1_tg['time'] = time


so4_a2 = xr.open_dataset(basepath+'ext/ar5_mam3_so4_a2_elev_1850_c090726.nc')
time = cal.MMDD2date(so4_a2.date)
contvolc_a2_tg = convert_molecules_to_tg(so4_a2,'contvolc', molecular_weight)
contvolc_a2_tg['time'] = time

so4_a1 = xr.open_dataset(basepath+'srf/ar5_mam3_so4_a1_surf_1850_c090726.nc')
time = cal.MMDD2date(so4_a1.date)
awb_tg = convert_molecules_to_tg(so4_a1,'emiss_awb', molecular_weight)
wst_tg = convert_molecules_to_tg(so4_a1,'emiss_wst', molecular_weight)
shp_tg = convert_molecules_to_tg(so4_a1,'emiss_shp', molecular_weight)
awb_tg['time'] = time
wst_tg['time'] = time
shp_tg['time'] = time


so4_a2 = xr.open_dataset(basepath+'srf/ar5_mam3_so4_a2_surf_1850_c090726.nc')
time = cal.MMDD2date(so4_a2.date)
dom_tg = convert_molecules_to_tg(so4_a2,'emiss_dom', molecular_weight)
tra_tg = convert_molecules_to_tg(so4_a2,'emiss_tra', molecular_weight)
dom_tg['time'] = time
tra_tg['time'] = time

#ene_tg = ene_tg.groupby('time.year').mean('time')
#ind_tg = ind_tg.groupby('time.year').mean('time')
#forestfire_tg = forestfire_tg.groupby('time.year').mean('time')
#grassfire_tg = grassfire_tg.groupby('time.year').mean('time')
#contvolc_a1_tg = contvolc_a1_tg.groupby('time.year').mean('time')
#contvolc_a2_tg = contvolc_a2_tg.groupby('time.year').mean('time')
#awb_tg = awb_tg.groupby('time.year').mean('time')
#wst_tg = wst_tg.groupby('time.year').mean('time')
#shp_tg = shp_tg.groupby('time.year').mean('time')
#dom_tg = dom_tg.groupby('time.year').mean('time')
#tra_tg = tra_tg.groupby('time.year').mean('time')

ene_tg = ene_tg.rename('ene_tg')
ind_tg = ind_tg.rename('ind_tg')
forestfire_tg = forestfire_tg.rename('forestfire_tg')
grassfire_tg = grassfire_tg.rename('grassfire_tg')
contvolc_a1_tg = contvolc_a1_tg.rename('contvolc_a1_tg')
contvolc_a2_tg = contvolc_a2_tg.rename('contvolc_a1_tg')
awb_tg = awb_tg.rename('awb_tg')
wst_tg = wst_tg.rename('wst_tg')
shp_tg = shp_tg.rename('shp_tg')
dom_tg = dom_tg.rename('dom_tg')
tra_tg = tra_tg.rename('tra_tg')


ene_tg.to_netcdf(pathout+'piControl_CESM1.nc')
ind_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
forestfire_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
grassfire_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
contvolc_a1_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
contvolc_a2_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
awb_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
wst_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
shp_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
dom_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')
tra_tg.to_netcdf(pathout+'piControl_CESM1.nc', mode='a')


#dat = xr.merge([ene_tg, ind_tg, forestfire_tg, grassfire_tg, contvolc_a1_tg,
#                contvolc_a2_tg, awb_tg, wst_tg, shp_tg, dom_tg, tra_tg])
#
#dat.to_netcdf(pathout+'piControl_CESM2.nc')












