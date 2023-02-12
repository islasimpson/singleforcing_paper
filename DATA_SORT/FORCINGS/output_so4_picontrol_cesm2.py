import xarray as xr
import numpy as np
import sys
import pandas as pd

from CASutils import calendar_utils as cal

avog = 6.022e23 # Avogadro's number
re = 6.3712e8 # Radius of the earth in cm

def convert_molecules_to_tg(dat,varname, molecular_weight=None):
    """ Convert surface emissions in molecules/cm2/s to Tg 

    """
    # Start with moleculre/cm2/s.  Convert from molecules to grams.  
    # Divide by Avogadro's number to convert from molecules to moles.  
    #Multiply by molcular weight in g/mol to end up with g/cm2/s
    if (molecular_weight):
        dat_g = molecular_weight*dat[varname]/avog
    else:
        dat_g = dat[varname].molecular_weight*dat[varname]/avog

    # convert from per s to per year
    dat_g_y = dat_g*365.*86400.

    if "altitude" in dat.dims:
        print('you have altitudes')
        dz = dat.altitude_int[1:dat.altitude_int.size].values - dat.altitude_int[0:dat.altitude_int.size-1].values
        dz = xr.DataArray(dz, coords=[dat.altitude], dims=['altitude'], name='dz')
        dz = dz*1000.*100. # convert to cm
        dat_g_y = (dat_g_y*dz).sum('altitude')
    if "lev" in dat[varname].dims:
        pressures = dat.hyam*1e5 + dat.hybm*dat.PS
        ipressures = dat.hyai*1e5 + dat.hybi*dat.PS
        dp = ipressures[1:ipressures.ilev.size,:,:,:] - np.array(ipressures[0:ipressures.ilev.size-1,:,:,:])
        dp = dp.rename({'ilev':'lev'})
        dp['lev'] = pressures.lev
        dz = (7.5/pressures)*dp # in km
        dat_g_y = (dat_g_y*dz).sum('lev')

#    if "altitude" in dat.dims:
#        print('you have altitudes')
#        dz = dat.altitude_int[1:dat.altitude_int.size].values - dat.altitude_int[0:dat.altitude_int.size-1].values
#        dz = xr.DataArray(dz, coords=[dat.altitude], dims=['altitude'], name='dz')
#        dz = dz*1000.*100. # convert to cm
#        dat_g_y = (dat_g_y*dz).sum('altitude')
#    else if  "lev" in dat[varname].dims:
#        pressures = dat.hyam*1e5 + dat.hybm*dat.PS
#        ipressures = dat.hyai*1e5 + dat.hybi*dat.PS
#        dp = ipressures[1:ipressures.ilev.size,:,:,:] - np.array(ipressures[0:ipressures.ilev.size-1,:,:,:])
#        dp = dp.rename({'ilev':'lev'})
#        dp['lev'] = pressures.lev
#        dz = (7.5/pressures)*dp # in km
#        dat_g_y = (dat_g_y*dz).sum(lev)


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

basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM2/piControl/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/"

so4_a1 = xr.open_dataset(basepath+'ext/emissions-cmip6_so4_a1_anthro-ag-ship_surface_1750-2015_0.9x1.25_c20170616.nc')
#time = cal.YYYYMMDD2date(so4_a1.date)
#so4_a1['time'] = time

ag_sol_was_tg = convert_molecules_to_tg(so4_a1,'emiss_ag_sol_was')
ag_sol_was_tg = ag_sol_was_tg.sel(time=slice("1850-01","1850-12"))
#ag_sol_was_tg = ag_sol_was_tg.groupby('time.year').mean('time')
ag_sol_was_tg = ag_sol_was_tg.rename('ag_sol_was_so4')
#ag_sol_was_tg = ag_sol_was_tg.sel(year=1850)


ship_tg = convert_molecules_to_tg(so4_a1,'emiss_shipping')
#ship_tg = ship_tg.groupby('time.year').mean('time')
ship_tg = ship_tg.sel(time=slice("1850-01","1850-12"))
ship_tg = ship_tg.rename('ship_so4')
#ship_tg = ship_tg.sel(year=1850)

so4_a1 = xr.open_dataset(basepath+
      'ext/emissions-cmip6_so4_a1_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc')
#time = cal.YYYYMMDD2date(so4_a1.date)
#so4_a1['time'] = time

contvolc_a1_tg = convert_molecules_to_tg(so4_a1,'contvolcano',molecular_weight=115.)
#contvolc_a1_tg = contvolc_a1_tg.groupby('time.year').mean('time')
contvolc_a1_tg = contvolc_a1_tg.sel(time=slice("1850-01","1850-12"))
contvolc_a1_tg = contvolc_a1_tg.rename('contvolc_a1_so4')
#contvolc_a1_tg = contvolc_a1_tg.sel(year=1850)

so4_a2 = xr.open_dataset(basepath+
      'ext/emissions-cmip6_so4_a2_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc')
#time = cal.YYYYMMDD2date(so4_a2.date)
#so4_a2['time'] = time

contvolc_a2_tg = convert_molecules_to_tg(so4_a2,'contvolcano', molecular_weight=115.)
#contvolc_a2_tg = contvolc_a2_tg.groupby('time.year').mean('time')
contvolc_a2_tg = contvolc_a2_tg.sel(time=slice("1850-01","1850-12"))
contvolc_a2_tg = contvolc_a2_tg.rename('contvolc_a2_so4')
#contvolc_a2_tg = contvolc_a2_tg.sel(year=1850)

so4_a1 = xr.open_dataset(basepath+
     'srf/emissions-cmip6_so4_a1_anthro-ag-ship_surface_1750-2015_0.9x1.25_c20170616.nc')
#time = cal.YYYYMMDD2date(so4_a1.date)
#so4_a1['time'] = time
ag_sol_was_srf_tg = convert_molecules_to_tg(so4_a1,'emiss_ag_sol_was')
#ag_sol_was_srf_tg = ag_sol_was_srf_tg.groupby('time.year').mean('time')
ag_sol_was_srf_tg = ag_sol_was_srf_tg.sel(time=slice("1850-01","1850-12"))
ag_sol_was_srf_tg = ag_sol_was_srf_tg.rename('ag_sol_was_so4_srf')

ship_srf_tg = convert_molecules_to_tg(so4_a1,'emiss_shipping')
#ship_srf_tg = ship_srf_tg.groupby('time.year').mean('time')
ship_srf_tg = ship_srf_tg.sel(time=slice("1850-01","1850-12"))
ship_srf_tg = ship_srf_tg.rename('ship_so4_srf')

so4_a1 = xr.open_dataset(basepath+
    'srf/emissions-cmip6_so4_a1_bb_surface_1750-2015_0.9x1.25_c20170322.nc')
#time = cal.YYYYMMDD2date(so4_a1.date)
#so4_a1['time'] = time
bb_tg = convert_molecules_to_tg(so4_a1, 'emiss_bb')
#bb_tg = bb_tg.groupby('time.year').mean('time')
bb_tg = bb_tg.sel(time=slice("1850-01","1850-12"))
bb_tg = bb_tg.rename('bb')


so4_a2 = xr.open_dataset(basepath+
    'srf/emissions-cmip6_so4_a2_anthro-res_surface_1750-2015_0.9x1.25_c20170616.nc')
#time = cal.YYYYMMDD2date(so4_a2.date)
#so4_a2['time'] = time
res_tran_tg = convert_molecules_to_tg(so4_a2,'emiss_res_tran')
#res_tran_tg = res_tran_tg.groupby('time.year').mean('time')
res_tran_tg = res_tran_tg.sel(time=slice("1850-01","1850-12"))
res_tran_tg = res_tran_tg.rename('res_tran')

time = pd.date_range("1850-01","1851-01",freq="M")
ag_sol_was_tg['time'] = time
ship_tg['time'] = time
contvolc_a1_tg['time'] = time
contvolc_a2_tg['time'] = time
ag_sol_was_srf_tg['time'] = time
ship_srf_tg['time'] = time
bb_tg['time'] = time
res_tran_tg['time'] = time


dat = xr.merge([ag_sol_was_tg,ship_tg,contvolc_a1_tg,contvolc_a2_tg,
                ag_sol_was_srf_tg,ship_srf_tg,bb_tg,res_tran_tg])

dat.to_netcdf(pathout+'piControl_CESM2.nc')



































