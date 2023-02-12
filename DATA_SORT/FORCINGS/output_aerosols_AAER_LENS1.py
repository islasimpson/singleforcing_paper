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

#----------------------------------------------------------------------------
#-----------------------------Black Carbon-----------------------------------
#----------------------------------------------------------------------------
anthro_bc = xr.open_dataset(basepath+'BC/srf/historical/ar5_mam3_bc_surf_1850-2005_c090804.nc')
molecular_weight = 12.
time = cal.YYYYMMDD2date(anthro_bc.date)
awb_tg = convert_molecules_to_tg(anthro_bc,'emiss_awb', molecular_weight)
dom_tg = convert_molecules_to_tg(anthro_bc,'emiss_dom', molecular_weight)
ene_tg = convert_molecules_to_tg(anthro_bc,'emiss_ene', molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_bc,'emiss_ind', molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_bc,'emiss_tra', molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_bc,'emiss_wst', molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_bc,'emiss_shp', molecular_weight)
anthro_bc_tg = awb_tg + dom_tg + ene_tg + ind_tg + tra_tg + wst_tg + shp_tg
anthro_bc_tg['time'] = time

anthro_bc_rcp = xr.open_dataset(basepath+'BC/srf/rcp85/RCP85_mam3_bc_surf_2000-2300_c20120214.nc')
anthro_bc_rcp = anthro_bc_rcp.where(anthro_bc_rcp.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_bc_rcp.date)
awb_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_awb', molecular_weight)
dom_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_dom', molecular_weight)
ene_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_ene', molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_ind', molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_tra', molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_wst', molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_bc_rcp,'emiss_shp', molecular_weight)
anthro_bc_tg_rcp = awb_tg + dom_tg + ene_tg + ind_tg + tra_tg + wst_tg + shp_tg
anthro_bc_tg_rcp['time'] = time

anthro_bc_tg = xr.concat([anthro_bc_tg,anthro_bc_tg_rcp], dim='time')
anthro_bc_tg = anthro_bc_tg.groupby('time.year').mean('time')
anthro_bc_tg = anthro_bc_tg.rename('anthro_bc')
anthro_bc_tg.attrs['units'] = 'Tg/y'

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#----------------------------------SO2---------------------------------------
#----------------------------------------------------------------------------
anthro_so2 = xr.open_dataset(basepath+'SO2/ext/historical/ar5_mam3_so2_elev_1850-2005_c090804.nc')
time = cal.YYYYMMDD2date(anthro_so2.date)
molecular_weight = 64.
ene_tg = convert_molecules_to_tg(anthro_so2,'emiss_ene', molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_so2,'emiss_ind', molecular_weight)

anthro_so2 = xr.open_dataset(basepath+'SO2/srf/historical/ar5_mam3_so2_surf_1850-2005_c090804.nc')
awb_tg = convert_molecules_to_tg(anthro_so2,'emiss_awb', molecular_weight)
dom_tg = convert_molecules_to_tg(anthro_so2,'emiss_dom', molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_so2,'emiss_tra', molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_so2,'emiss_wst', molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_so2,'emiss_shp', molecular_weight)

anthro_so2_tg = ene_tg + ind_tg + awb_tg + dom_tg + tra_tg + wst_tg + shp_tg
anthro_so2_tg['time'] = time

anthro_so2 = xr.open_dataset(basepath+'SO2/ext/rcp85/RCP85_mam3_so2_elev_2000-2300_c20120214.nc')
anthro_so2 = anthro_so2.where(anthro_so2.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_so2.date)
molecular_weight = 64.
ene_tg = convert_molecules_to_tg(anthro_so2,'emiss_ene', molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_so2,'emiss_ind', molecular_weight)

anthro_so2 = xr.open_dataset(basepath+'SO2/srf/rcp85/RCP85_mam3_so2_surf_2000-2300_c20120214.nc')
anthro_so2 = anthro_so2.where(anthro_so2.date < 21010000, drop=True)
awb_tg = convert_molecules_to_tg(anthro_so2,'emiss_awb', molecular_weight)
dom_tg = convert_molecules_to_tg(anthro_so2,'emiss_dom', molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_so2,'emiss_tra', molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_so2,'emiss_wst', molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_so2,'emiss_shp', molecular_weight)

anthro_so2_tg_rcp = ene_tg + ind_tg + awb_tg + dom_tg + tra_tg + wst_tg + shp_tg
anthro_so2_tg_rcp['time'] = time

anthro_so2_tg = xr.concat([anthro_so2_tg,anthro_so2_tg_rcp], dim='time')

anthro_so2_tg = anthro_so2_tg.groupby('time.year').mean('time')
anthro_so2_tg = anthro_so2_tg.rename('anthro_so2')
anthro_so2_tg.attrs['units'] = 'Tg/y'
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#----------------------------------SO4----------------------------------------
#-----------------------------------------------------------------------------
anthro_so4_a1 = xr.open_dataset(basepath+'SO4/ext/historical/ar5_mam3_so4_a1_elev_1850-2005_c090804.nc')
time = cal.YYYYMMDD2date(anthro_so4_a1.date)
molecular_weight = 115.
ene_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_ene',molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_ind',molecular_weight)

anthro_so4_a1 = xr.open_dataset(basepath+'SO4/srf/historical/ar5_mam3_so4_a1_surf_1850-2005_c090804.nc')
awb_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_awb',molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_wst',molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_shp',molecular_weight)

anthro_so4_a2 = xr.open_dataset(basepath+'SO4/srf/historical/ar5_mam3_so4_a2_surf_1850-2005_c090804.nc')
dom_tg = convert_molecules_to_tg(anthro_so4_a2,'emiss_dom',molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_so4_a2,'emiss_tra', molecular_weight)

anthro_so4_tg = ene_tg + ind_tg + awb_tg + wst_tg + shp_tg + dom_tg + tra_tg
anthro_so4_tg['time'] = time

anthro_so4_a1 = xr.open_dataset(basepath+'SO4/ext/rcp85/RCP85_mam3_so4_a1_elev_2000-2300_c20120214.nc')
anthro_so4_a1 = anthro_so4_a1.where(anthro_so4_a1.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_so4_a1.date)
molecular_weight = 115.
ene_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_ene',molecular_weight)
ind_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_ind',molecular_weight)

anthro_so4_a1 = xr.open_dataset(basepath+'SO4/srf/rcp85/RCP85_mam3_so4_a1_surf_2000-2300_c20120214.nc')
anthro_so4_a1 = anthro_so4_a1.where(anthro_so4_a1.date < 21010000, drop=True)
awb_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_awb',molecular_weight)
wst_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_wst',molecular_weight)
shp_tg = convert_molecules_to_tg(anthro_so4_a1,'emiss_shp',molecular_weight)

anthro_so4_a2 = xr.open_dataset(basepath+'SO4/srf/rcp85/RCP85_mam3_so4_a2_surf_2000-2300_c20120214.nc')
anthro_so4_a2 = anthro_so4_a2.where(anthro_so4_a2.date < 21010000, drop=True)
dom_tg = convert_molecules_to_tg(anthro_so4_a2,'emiss_dom',molecular_weight)
tra_tg = convert_molecules_to_tg(anthro_so4_a2,'emiss_tra', molecular_weight)

anthro_so4_tg_rcp = ene_tg + ind_tg + awb_tg + wst_tg + shp_tg + dom_tg + tra_tg
anthro_so4_tg_rcp['time'] = time

anthro_so4_tg = xr.concat([anthro_so4_tg, anthro_so4_tg_rcp], dim='time')

anthro_so4_tg = anthro_so4_tg.groupby('time.year').mean('time')
anthro_so4_tg = anthro_so4_tg.rename('anthro_so4')
anthro_so4_tg.attrs['units'] = 'Tg/y'


#-----------------------------------------------------------------------------
#--------------------------------------SOAG-----------------------------------
#-----------------------------------------------------------------------------
anthro_soag = xr.open_dataset(basepath+'SOAG/srf/historical/'+\
         'ar5_mam3_soag_1.5_surf_1850-2005_c100429.nc')
anthro_soag = anthro_soag.where(anthro_soag.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_soag.date)
molecular_weight=12.
anthro_soag_bigalk = convert_molecules_to_tg(0.05*anthro_soag,'SOAG_BIGALK',molecular_weight)
#anthro_soag_bigalk = 0.05*anthro_soag_bigalk # multiplying by the yield
anthro_soag_bigene = convert_molecules_to_tg(0.05*anthro_soag,'SOAG_BIGENE',molecular_weight)
#anthro_soag_bigene = 0.05*anthro_soag_bigene
anthro_soag_terpene = convert_molecules_to_tg(0.25*anthro_soag,'SOAG_TERPENE',molecular_weight)
#anthro_soag_terpene = 0.25*anthro_soag_terpene
anthro_soag_toluene = convert_molecules_to_tg(0.15*anthro_soag,'SOAG_TOLUENE',molecular_weight)
#anthro_soag_toluene = 0.15*anthro_soag_toluene
anthro_soag_tg = anthro_soag_bigalk + anthro_soag_bigene + \
               anthro_soag_terpene + anthro_soag_toluene
#anthro_soag_tg = anthro_soag_bigalk + anthro_soag_bigene + anthro_soag_toluene
anthro_soag_tg['time'] = time

anthro_soag = xr.open_dataset(basepath+'SOAG/srf/rcp85/'+\
         'RCP85_soag_1.5_surf_2000-2300_c20120214.nc')
anthro_soag = anthro_soag.where(anthro_soag.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_soag.date)
molecular_weight=12.
anthro_soag_bigalk = convert_molecules_to_tg(0.05*anthro_soag,'SOAG_BIGALK',molecular_weight)
#anthro_soag_bigalk = 0.05*anthro_soag_bigalk
anthro_soag_bigene = convert_molecules_to_tg(0.05*anthro_soag,'SOAG_BIGENE',molecular_weight)
#anthro_soag_bigene = 0.05*anthro_soag_bigene
anthro_soag_terpene = convert_molecules_to_tg(0.25*anthro_soag,'SOAG_TERPENE',molecular_weight)
#anthro_soag_terpene = 0.25*anthro_soag_terpene
anthro_soag_toluene = convert_molecules_to_tg(0.15*anthro_soag,'SOAG_TOLUENE',molecular_weight)
#anthro_soag_toluene = 0.15*anthro_soag_toluene
anthro_soag_tg_rcp85 = anthro_soag_bigalk + anthro_soag_bigene + \
               anthro_soag_terpene + anthro_soag_toluene
#anthro_soag_tg_rcp85 = anthro_soag_bigalk + anthro_soag_bigene + anthro_soag_toluene

anthro_soag_tg_rcp85['time'] = time

anthro_soag_tg = xr.concat([anthro_soag_tg, anthro_soag_tg_rcp85], dim='time')
anthro_soag_tg = anthro_soag_tg.groupby('time.year').mean('time')
anthro_soag_tg = anthro_soag_tg.rename('anthro_soag')
anthro_soag_tg.attrs['units'] = 'Tg/y'


#-------------------------------------------------------------------------------
#------------------------------------POM----------------------------------------
#-------------------------------------------------------------------------------
anthro_pom = xr.open_dataset(basepath+'pom/srf/historical/'+\
          'ar5_mam3_oc_surf_1850-2005_c090804.nc')
anthro_pom = anthro_pom.where(anthro_pom.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_pom.date)
molecular_weight=12.
anthro_pom_awb = convert_molecules_to_tg(anthro_pom,'emiss_awb',molecular_weight)
anthro_pom_dom = convert_molecules_to_tg(anthro_pom,'emiss_dom',molecular_weight)
anthro_pom_ene = convert_molecules_to_tg(anthro_pom,'emiss_ene',molecular_weight)
anthro_pom_ind = convert_molecules_to_tg(anthro_pom,'emiss_ind',molecular_weight)
anthro_pom_tra = convert_molecules_to_tg(anthro_pom,'emiss_tra',molecular_weight)
anthro_pom_shp = convert_molecules_to_tg(anthro_pom,'emiss_shp',molecular_weight)

anthro_pom_tg = anthro_pom_awb + anthro_pom_dom + anthro_pom_ene + \
                    anthro_pom_ind + anthro_pom_tra + anthro_pom_shp
anthro_pom_tg['time'] = time


anthro_pom = xr.open_dataset(basepath+'pom/srf/rcp85/'+\
        'RCP85_mam3_oc_surf_2000-2300_c20120214.nc')
anthro_pom = anthro_pom.where(anthro_pom.date < 21010000, drop=True)
time = cal.YYYYMMDD2date(anthro_pom.date)
molecular_weight=12.
anthro_pom_awb = convert_molecules_to_tg(anthro_pom,'emiss_awb',molecular_weight)
anthro_pom_dom = convert_molecules_to_tg(anthro_pom,'emiss_dom',molecular_weight)
anthro_pom_ene = convert_molecules_to_tg(anthro_pom,'emiss_ene',molecular_weight)
anthro_pom_ind = convert_molecules_to_tg(anthro_pom,'emiss_ind',molecular_weight)
anthro_pom_tra = convert_molecules_to_tg(anthro_pom,'emiss_tra',molecular_weight)
anthro_pom_shp = convert_molecules_to_tg(anthro_pom,'emiss_shp',molecular_weight)

anthro_pom_tg_rcp85 = anthro_pom_awb + anthro_pom_dom + anthro_pom_ene + \
                      anthro_pom_ind + anthro_pom_tra + anthro_pom_shp
anthro_pom_tg_rcp85['time'] = time
anthro_pom_tg = xr.concat([anthro_pom_tg, anthro_pom_tg_rcp85], dim='time')
anthro_pom_tg = anthro_pom_tg.groupby('time.year').mean('time')
anthro_pom_tg = anthro_pom_tg.rename('anthro_pom')
anthro_pom_tg.attrs['units'] = 'Tg/y'

# interpolate onto a common year axis.
yuse = np.arange(1850,2050,1)
anthro_bc_tg = anthro_bc_tg.interp(year=yuse)
anthro_so2_tg = anthro_so2_tg.interp(year=yuse)
anthro_so4_tg = anthro_so4_tg.interp(year=yuse)
anthro_soag_tg = anthro_soag_tg.interp(year=yuse)
anthro_pom_tg = anthro_pom_tg.interp(year=yuse)


#---------------------------------------------------------
anthro_bc_tg.to_netcdf(pathout+'AAER_LENS1.nc')
anthro_so2_tg.to_netcdf(pathout+'AAER_LENS1.nc', mode='a')
anthro_so4_tg.to_netcdf(pathout+'AAER_LENS1.nc', mode='a')
anthro_soag_tg.to_netcdf(pathout+'AAER_LENS1.nc', mode='a')
anthro_pom_tg.to_netcdf(pathout+'AAER_LENS1.nc', mode='a')















