import xarray as xr
import numpy as np
import sys

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

basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM2/BMB/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/"

#---------------------------------------------------------------------------------
#-----------------------------Black Carbon----------------------------------------
#---------------------------------------------------------------------------------
bmb_bc = xr.open_dataset(basepath+'BC/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_bc_a4_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_bc_tg = convert_molecules_to_tg(bmb_bc,'emiss_bb')
bmb_bc_tg = bmb_bc_tg.groupby('time.year').mean('time')
bmb_bc_tg = bmb_bc_tg.rename('bmb_bc')
bmb_bc_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
#--------------------------------DMS----------------------------------------------
#---------------------------------------------------------------------------------
bmb_dms = xr.open_dataset(basepath+'DMS/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_DMS_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_dms_tg = convert_molecules_to_tg(bmb_dms,'emiss_bb')
bmb_dms_tg = bmb_dms_tg.groupby('time.year').mean('time')
bmb_dms_tg = bmb_dms_tg.rename('bmb_dms')
bmb_dms_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#---------------------------------pom---------------------------------------------
#---------------------------------------------------------------------------------
bmb_pom = xr.open_dataset(basepath+'/pom/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_pom_a4_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_pom_tg = convert_molecules_to_tg(bmb_pom,'emiss_bb')
bmb_pom_tg = bmb_pom_tg.groupby('time.year').mean('time')
bmb_pom_tg = bmb_pom_tg.rename('bmb_pom')
bmb_pom_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
#---------------------------------SO2---------------------------------------------
#---------------------------------------------------------------------------------
bmb_so2 = xr.open_dataset(basepath+'/SO2/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_SO2_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_so2_tg = convert_molecules_to_tg(bmb_so2,'emiss_bb')
bmb_so2_tg = bmb_so2_tg.groupby('time.year').mean('time')
bmb_so2_tg = bmb_so2_tg.rename('bmb_so2')
bmb_so2_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
#---------------------------------SO4---------------------------------------------
#---------------------------------------------------------------------------------
bmb_so4 = xr.open_dataset(basepath+'/SO4/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_so4_a1_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_so4_tg = convert_molecules_to_tg(bmb_so4,'emiss_bb')
bmb_so4_tg = bmb_so4_tg.groupby('time.year').mean('time')
bmb_so4_tg = bmb_so4_tg.rename('bmb_so4')
bmb_so4_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
#---------------------------------SOAG---------------------------------------------
#---------------------------------------------------------------------------------
bmb_soag = xr.open_dataset(basepath+'/SOAG/srf/historical/'+\
'emissions-cmip6-ScenarioMIP_IAMC-AIM-ssp370-1-1_smoothed_SOAGx1.5_bb_surface_mol_175001-210101_0.9x1.25_c20201016.nc')
bmb_soag_tg = convert_molecules_to_tg(bmb_soag,'emiss_bb')
bmb_soag_tg = bmb_soag_tg.groupby('time.year').mean('time')
bmb_soag_tg = bmb_soag_tg.rename('bmb_soag')
bmb_soag_tg.attrs['units'] = 'Tg/y'
#---------------------------------------------------------------------------------


#interpolate onto a common year axis
yuse = np.arange(1850,2050,1)
bmb_bc_tg = bmb_bc_tg.interp(year=yuse)
bmb_dms_tg = bmb_dms_tg.interp(year=yuse)
bmb_pom_tg = bmb_pom_tg.interp(year=yuse)
bmb_so2_tg = bmb_so2_tg.interp(year=yuse)
bmb_so4_tg = bmb_so4_tg.interp(year=yuse)
bmb_soag_tg = bmb_soag_tg.interp(year=yuse)

bmb_bc_tg.to_netcdf(pathout+'BMB.nc')
bmb_dms_tg.to_netcdf(pathout+'BMB.nc', mode='a')
bmb_pom_tg.to_netcdf(pathout+'BMB.nc', mode='a')
bmb_so2_tg.to_netcdf(pathout+'BMB.nc', mode='a')
bmb_so4_tg.to_netcdf(pathout+'BMB.nc', mode='a')
bmb_soag_tg.to_netcdf(pathout+'BMB.nc', mode='a')

