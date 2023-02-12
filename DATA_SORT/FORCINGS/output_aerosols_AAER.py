import xarray as xr
import numpy as np
import sys

from CASutils import calendar_utils as cal

avog = 6.022e23 # Avogadro's number
re = 6.3712e8 # Radius of the earth in cm

def convert_molecules_to_tg(dat,varname):
    """ Convert surface emissions in molecules/cm2/s to Tg 

    """
    # Start with moleculre/cm2/s.  Convert from molecules to grams.  
    # Divide by Avogadro's number to convert from molecules to moles.  
    #Multiply by molcular weight in g/mol to end up with g/cm2/s
    dat_g = dat[varname].molecular_weight*dat[varname]/avog

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
    
basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM2/AAER/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/"

#---------------------------------------------------------------------
#---------------------------Black carbon------------------------------
#---------------------------------------------------------------------
anthro_bc = xr.open_dataset(basepath+'BC/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_bc_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc')
#---building the time axis from date as can't trust the time variable in these files.
time = cal.YYYYMMDD2date(anthro_bc.date)
anthro_bc['time'] = time


#----Anthropogenic emissions (time evolving)
anthro_bc_tg = convert_molecules_to_tg(anthro_bc,'emiss_anthro')
anthro_bc_tg = anthro_bc_tg.groupby('time.year').mean('time')
anthro_bc_tg = anthro_bc_tg.rename('anthro_bc')
anthro_bc_tg.attrs['units']='Tg/y'
#----------------------------------------------------------------------




#----------------------------------------------------------------------
#----------------------------------SO2---------------------------------
#----------------------------------------------------------------------
eandi_so2 = xr.open_dataset(basepath+'SO2/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SO2_anthro-ene_surface_mol_175001-210101_0.9x1.25_c20190222.nc')
time = cal.YYYYMMDD2date(eandi_so2.date)
eandi_so2['time'] = time

agships_so2 = xr.open_dataset(basepath+'SO2/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SO2_anthro-ag-ship-res_surface_mol_175001-210101_0.9x1.25_c20200924.nc')
time = cal.YYYYMMDD2date(agships_so2.date)
agships_so2['time'] = time


#----Energy and Industrial emissions (time evolving)
eandi_so2_tg = convert_molecules_to_tg(eandi_so2,'emiss_ene_ind')
eandi_so2_tg = eandi_so2_tg.groupby('time.year').mean('time')
eandi_so2_tg = eandi_so2_tg.rename('ene_ind_so2')
eandi_so2_tg.attrs['units'] = 'Tg/y'

#---Agriculture, solvents and waste (time evolving)
ag_so2_tg = convert_molecules_to_tg(agships_so2,'emiss_ag_sol_was')
ag_so2_tg = ag_so2_tg.groupby('time.year').mean('time')
ag_so2_tg = ag_so2_tg.rename('ag_so2')
ag_so2_tg.attrs['units'] = 'Tg/y'

#---Residential and transport emissions (time evolving)
restran_so2_tg = convert_molecules_to_tg(agships_so2,'emiss_res_tran')
restran_so2_tg = restran_so2_tg.groupby('time.year').mean('time')
restran_so2_tg = restran_so2_tg.rename('restran_so2')
restran_so2_tg.attrs['units'] = 'Tg/y'

#---Ships emissions (time evolving)
ship_so2_tg = convert_molecules_to_tg(agships_so2,'emiss_ship')
ship_so2_tg = ship_so2_tg.groupby('time.year').mean('time')
ship_so2_tg = ship_so2_tg.rename('ship_so2')
ship_so2_tg.attrs['units'] = 'Tg/y'
#--------------------------------------------------------

#--------------------------------------------------------
#-----------------------SO4------------------------------
#--------------------------------------------------------
eneind_so4_hist = xr.open_dataset(basepath+'SO4/ext/historical/emissions-cmip6_so4_a1_anthro-ene_vertical_1750-2015_0.9x1.25_c20170616.nc')
time = cal.YYYYMMDD2date(eneind_so4_hist.date)
eneind_so4_hist['time'] = time
eneind_so4_hist = eneind_so4_hist.sel(time=slice("1750-01-01","2014-12-31"))

eneind_so4_ssp = xr.open_dataset(basepath+'SO4/ext/SSP370/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a1_anthro-ene_vertical_mol_175001-210101_0.9x1.25_c20190222.nc')
time = cal.YYYYMMDD2date(eneind_so4_ssp.date)
eneind_so4_ssp['time'] = time

eneind_so4_ssp = eneind_so4_ssp.sel(time=slice("2015-01-01","2100-12-31"))
eneind_so4 = xr.concat([eneind_so4_hist, eneind_so4_ssp], dim='time')

agships_so4 = xr.open_dataset(basepath+'SO4/srf/historical/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a1_anthro-ag-ship_surface_mol_175001-210101_0.9x1.25_c20200924.nc')
time = cal.YYYYMMDD2date(agships_so4.date)
agships_so4['time'] = time

restran_so4 = xr.open_dataset(basepath+'SO4/srf/historical/emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_so4_a2_anthro-res_surface_mol_175001-210101_0.9x1.25_c20200924.nc')
time = cal.YYYYMMDD2date(restran_so4.date)
restran_so4['time'] = time


#----Energy and Industrial emissions (time evolving)
eandi_so4_tg = convert_molecules_to_tg(eneind_so4,'emiss_ene_ind')
eandi_so4_tg = eandi_so4_tg.groupby('time.year').mean('time')
eandi_so4_tg = eandi_so4_tg.rename('ene_ind_so4')
eandi_so4_tg.attrs['units'] = 'Tg/y'

#----Agriculture, solvents and waste (time evolving)
ag_so4_tg = convert_molecules_to_tg(agships_so4,'emiss_ag_sol_was')
ag_so4_tg = ag_so4_tg.groupby('time.year').mean('time')
ag_so4_tg = ag_so4_tg.rename('ag_so4')
ag_so4_tg.attrs['units'] = 'Tg/y'

#----Residential and transport emissions (time evolving)
restran_so4_tg = convert_molecules_to_tg(restran_so4,'emiss_res_tran')
restran_so4_tg = restran_so4_tg.groupby('time.year').mean('time')
restran_so4_tg = restran_so4_tg.rename('restran_so4')
restran_so4_tg.attrs['units'] = 'Tg/y'

#---Ships emissions (time evolving)
ship_so4_tg = convert_molecules_to_tg(agships_so4,'emiss_shipping')
ship_so4_tg = ship_so4_tg.groupby('time.year').mean('time')
ship_so4_tg = ship_so4_tg.rename('ship_so4')
ship_so4_tg.attrs['units'] = 'Tg/y'
#--------------------------------------------------------

#--------------------------------------------------------
#----------------------POM-------------------------------
#--------------------------------------------------------
pom = xr.open_dataset(basepath+"/pom/srf/historical/"+\
"emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_pom_a4_anthro_surface_mol_175001-210101_0.9x1.25_c20190222.nc")
time = cal.YYYYMMDD2date(pom.date)
pom['time'] = time


pom_tg = convert_molecules_to_tg(pom,'emiss_anthro')
pom_tg = pom_tg.groupby('time.year').mean('time')
pom_tg = pom_tg.rename('pom')
pom_tg.attrs['units'] = 'Tg/y'
#--------------------------------------------------------

#--------------------------------------------------------
#---------------------SOAG-------------------------------
#--------------------------------------------------------
soag = xr.open_dataset(basepath+"/SOAG/srf/historical/"+\
"emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_SOAGx1.5_anthro_surface_mol_175001-210101_0.9x1.25_c20200403.nc")
time = cal.YYYYMMDD2date(soag.date)
soag['time'] = time

soag_tg = convert_molecules_to_tg(soag,'emiss_anthro')
soag_tg = soag_tg.groupby('time.year').mean('time')
soag_tg = soag_tg.rename('soag')
soag_tg.attrs['units'] = 'Tg/y'
#-------------------------------------------------------

# interpolate onto a common year axis.  Different years are present for
# different forcers into the future
yuse = np.arange(1850,2050,1) 
anthro_bc_tg =anthro_bc_tg.interp(year=yuse)
eandi_so2_tg = eandi_so2_tg.interp(year=yuse)
ag_so2_tg = ag_so2_tg.interp(year=yuse)
restran_so2_tg = restran_so2_tg.interp(year=yuse)
ship_so2_tg = ship_so2_tg.interp(year=yuse)
eandi_so4_tg = eandi_so4_tg.interp(year=yuse)
ag_so4_tg = ag_so4_tg.interp(year=yuse)
restran_so4_tg = restran_so4_tg.interp(year=yuse)
ship_so4_tg = ship_so4_tg.interp(year=yuse)
pom_tg = pom_tg.interp(year=yuse)
soag_tg = soag_tg.interp(year=yuse)

anthro_bc_tg.to_netcdf(pathout+'AAER.nc')
eandi_so2_tg.to_netcdf(pathout+'AAER.nc', mode='a')
ag_so2_tg.to_netcdf(pathout+'AAER.nc', mode='a')
restran_so2_tg.to_netcdf(pathout+'AAER.nc', mode='a')
ship_so2_tg.to_netcdf(pathout+'AAER.nc', mode='a')
eandi_so4_tg.to_netcdf(pathout+'AAER.nc', mode='a')
ag_so4_tg.to_netcdf(pathout+'AAER.nc', mode='a')
restran_so4_tg.to_netcdf(pathout+'AAER.nc', mode='a')
ship_so4_tg.to_netcdf(pathout+'AAER.nc', mode='a')
pom_tg.to_netcdf(pathout+'AAER.nc', mode='a')
soag_tg.to_netcdf(pathout+'AAER.nc', mode='a')


