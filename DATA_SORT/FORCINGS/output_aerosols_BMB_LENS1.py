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


basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM1/AAER/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/"


#------------------------------------------------------------------------------
#---------------------------Black Carbon---------------------------------------
#------------------------------------------------------------------------------
bmb_bc = xr.open_dataset(basepath+'BC/ext/historical/'+\
    'ar5_mam3_bc_elev_1850-2005_c090804.nc')
bmb_bc = bmb_bc.where(bmb_bc.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_bc.date)
molecular_weight=12.
bmb_bc_forest_tg = convert_molecules_to_tg(bmb_bc,'forestfire', molecular_weight)
bmb_bc_grass_tg = convert_molecules_to_tg(bmb_bc,'grassfire', molecular_weight)
bmb_bc_tg = bmb_bc_forest_tg + bmb_bc_grass_tg
bmb_bc_tg['time'] = time

bmb_bc = xr.open_dataset(basepath+'BC/ext/rcp85/'+\
    'RCP85_mam3_bc_elev_2000-2300_c20120214.nc')
bmb_bc = bmb_bc.where(bmb_bc.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_bc.date)
bmb_bc_forest_tg = convert_molecules_to_tg(bmb_bc, 'forestfire', molecular_weight)
bmb_bc_grass_tg = convert_molecules_to_tg(bmb_bc, 'grassfire', molecular_weight)
bmb_bc_tg_rcp85 = bmb_bc_forest_tg + bmb_bc_grass_tg
bmb_bc_tg_rcp85['time'] = time

bmb_bc_tg = xr.concat([bmb_bc_tg, bmb_bc_tg_rcp85], dim='time')
bmb_bc_tg = bmb_bc_tg.groupby('time.year').mean('time')
bmb_bc_tg = bmb_bc_tg.rename('bmb_bc')
bmb_bc_tg.attrs['units'] = 'Tg/y'

bmb_bc_tg.to_netcdf(pathout+'BMB_LENS1.nc')

#---------------------------------------------------------------------------------
#----------------------------POM--------------------------------------------------
#---------------------------------------------------------------------------------
bmb_pom = xr.open_dataset(basepath+'pom/ext/historical/'+\
      'ar5_mam3_oc_elev_1850-2005_c090804.nc')
bmb_pom = bmb_pom.where(bmb_pom.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_pom.date)
molecular_weight=12.
bmb_pom_forest_tg = convert_molecules_to_tg(bmb_pom,'forestfire',molecular_weight)
bmb_pom_grass_tg = convert_molecules_to_tg(bmb_pom, 'grassfire',molecular_weight)
bmb_pom_tg = bmb_pom_forest_tg + bmb_pom_grass_tg
bmb_pom_tg['time'] = time

bmb_pom = xr.open_dataset(basepath+'pom/ext/rcp85/'+\
     'RCP85_mam3_oc_elev_2000-2300_c20120214.nc')
bmb_pom = bmb_pom.where(bmb_pom.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_pom.date)
bmb_pom_forest_tg = convert_molecules_to_tg(bmb_pom,'forestfire',molecular_weight)
bmb_pom_grass_tg = convert_molecules_to_tg(bmb_pom,'grassfire',molecular_weight)
bmb_pom_tg_rcp85 = bmb_pom_forest_tg + bmb_pom_grass_tg
bmb_pom_tg_rcp85['time'] = time

bmb_pom_tg = xr.concat([bmb_pom_tg, bmb_pom_tg_rcp85], dim='time')
bmb_pom_tg = bmb_pom_tg.groupby('time.year').mean('time')
bmb_pom_tg = bmb_pom_tg.rename('bmb_pom')
bmb_pom_tg.attrs['units'] = 'Tg/y'

bmb_pom_tg.to_netcdf(pathout+'BMB_LENS1.nc', mode='a')






#---------------------------------------------------------------------------------
#-----------------------------SO4-------------------------------------------------
#---------------------------------------------------------------------------------

bmb_so4_a1 = xr.open_dataset(basepath+'SO4/ext/historical/'+\
   'ar5_mam3_so4_a1_elev_1850-2005_c090804.nc')
bmb_so4_a1 = bmb_so4_a1.where( bmb_so4_a1.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_so4_a1.date)
molecular_weight=115.

bmb_so4_a1_forest_tg = convert_molecules_to_tg(bmb_so4_a1, 'forestfire', molecular_weight)
bmb_so4_a1_grass_tg = convert_molecules_to_tg(bmb_so4_a1, 'grassfire', molecular_weight)
bmb_so4_a1_tg = bmb_so4_a1_forest_tg + bmb_so4_a1_grass_tg
bmb_so4_a1_tg['time'] = time

bmb_so4_hist = bmb_so4_a1_tg 


bmb_so4_a1_rcp85 = xr.open_dataset(basepath+'SO4/ext/rcp85/'+\
   'RCP85_mam3_so4_a1_elev_2000-2300_c20120214.nc')
bmb_so4_a1_rcp85 = bmb_so4_a1_rcp85.where( bmb_so4_a1_rcp85.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_so4_a1_rcp85.date)
bmb_so4_a1_forest_tg_rcp85 = convert_molecules_to_tg(bmb_so4_a1_rcp85, 'forestfire', molecular_weight)
bmb_so4_a1_grass_tg_rcp85 = convert_molecules_to_tg(bmb_so4_a1_rcp85, 'grassfire', molecular_weight)
bmb_so4_a1_tg_rcp85 = bmb_so4_a1_forest_tg_rcp85 + bmb_so4_a1_grass_tg_rcp85
bmb_so4_a1_tg_rcp85['time'] = time

bmb_so4_rcp85 = bmb_so4_a1_tg_rcp85 

bmb_so4_tg = xr.concat([bmb_so4_hist, bmb_so4_rcp85], dim='time')
bmb_so4_tg = bmb_so4_tg.groupby('time.year').mean('time')
bmb_so4_tg = bmb_so4_tg.rename('bmb_so4')
bmb_so4_tg.attrs['units'] = 'Tg/y'

bmb_so4_tg.to_netcdf(pathout+'BMB_LENS1.nc', mode='a')


#-----------------------------------------------------------------------------------
#--------------------------------SO2------------------------------------------------
#-----------------------------------------------------------------------------------
bmb_so2 = xr.open_dataset(basepath+'SO2/ext/historical/'+\
   'ar5_mam3_so2_elev_1850-2005_c090804.nc')
bmb_so2 = bmb_so2.where( bmb_so2.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_so2.date)
molecular_weight=64.

bmb_so2_forest_tg = convert_molecules_to_tg(bmb_so2, 'forestfire', molecular_weight)
bmb_so2_grass_tg = convert_molecules_to_tg(bmb_so2, 'grassfire', molecular_weight)
bmb_so2_tg = bmb_so2_forest_tg + bmb_so2_grass_tg
bmb_so2_tg['time'] = time

bmb_so2_hist = bmb_so2_tg


bmb_so2_rcp85 = xr.open_dataset(basepath+'SO2/ext/rcp85/'+\
   'RCP85_mam3_so2_elev_2000-2300_c20120214.nc')
bmb_so2_rcp85 = bmb_so2_rcp85.where( bmb_so2_rcp85.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(bmb_so2_rcp85.date)
bmb_so2_forest_tg_rcp85 = convert_molecules_to_tg(bmb_so2_rcp85, 'forestfire', molecular_weight)
bmb_so2_grass_tg_rcp85 = convert_molecules_to_tg(bmb_so2_rcp85, 'grassfire', molecular_weight)
bmb_so2_tg_rcp85 = bmb_so2_forest_tg_rcp85 + bmb_so2_grass_tg_rcp85
bmb_so2_tg_rcp85['time'] = time

bmb_so2_rcp85 = bmb_so2_tg_rcp85 

bmb_so2_tg = xr.concat([bmb_so2_hist, bmb_so2_rcp85], dim='time')
bmb_so2_tg = bmb_so2_tg.groupby('time.year').mean('time')
bmb_so2_tg = bmb_so2_tg.rename('bmb_so2')
bmb_so2_tg.attrs['units'] = 'Tg/y'

bmb_so2_tg.to_netcdf(pathout+'BMB_LENS1.nc', mode='a')





























