# Output annual means for the second 50 members of the CESM2 large ensemble
import xarray as xr
import numpy as np
import glob
import sys
import dask

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

import importlib
importlib.reload(lens)

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"
topdir="/project/mojave/cesm2/LENS/atm/month_1/"

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def calculateweights(darray):
    """Calculate area weights"""
    darray = landfrac
    a = 6.371e6
    latrad = np.deg2rad(darray.lat)
    dphi = xr.DataArray(np.zeros([latrad.size]), coords=[darray.lat], dims=['lat'], name='dphi')
    dphimid = [ ((latrad.isel(lat=i+1).values-latrad.isel(lat=i).values)/2. +  (latrad.isel(lat=i).values - latrad.isel(lat=i-1).values)/2.) 
            for i in np.arange(1,latrad.size-1,1)  ]
    dphi[1:darray.lat.size-1] = dphimid
    dphi[0] = (latrad.isel(lat=1).values - latrad.isel(lat=0).values)/2. + (latrad.isel(lat=0).values-(np.deg2rad(-90.)))
    dphi[darray.lat.size-1] = (latrad.isel(lat=latrad.size-1).values - latrad.isel(lat=latrad.size-2).values)/2. + \
                                   np.deg2rad(90) - latrad.isel(lat=latrad.size-1).values
    lonrad = np.deg2rad(darray.lon)
    dlon = lonrad[1].values - lonrad[0].values
    
    cosphi = xr.DataArray(np.zeros([latrad.size]), coords=[darray.lat], dims=['lat'], name='cosphi')
    cosphi[:] = np.cos(np.deg2rad(darray.lat.values))
        
    weights = xr.DataArray(np.ones([latrad.size, lonrad.size]), coords=[darray.lat, darray.lon], dims=['lat','lon'], name='weights')
    weights = a**2 * weights*cosphi*dphi*dlon
    return weights


#varnames=['AODVIS','FLNT','FSNT','FSNS','FSDS']
varnames=['FSNS','FSDS']

memstr = lens.lens2memnamegen_second50(50)

for varname in varnames:
    print(varname)
    filelist=[sorted(glob.glob(topdir+varname+"/*.BHISTsmbb.*"+imem+"*.nc"))+
              sorted(glob.glob(topdir+varname+"/*.BSSP370smbb.*"+imem+"*.nc"))
              for imem in memstr ]
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[[varname,'time_bnds']] for i in filelist]

    dat = xr.concat(members, dim='M', join='override', coords='minimal')
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    dat['lon'] = landfrac.lon ; dat['lat'] = landfrac.lat
   
    weights = calculateweights(dat)
 
    landfrac_stack = landfrac.landfrac.stack(z=("lon","lat"))
    weights_stack = weights.stack(z=("lon","lat"))
    dat_stack = dat.stack(z=("lon","lat"))

    dat_land = dat_stack.where(landfrac_stack > 0.5, drop=True)
    weights_land = weights_stack.where(landfrac_stack > 0.5, drop=True)
    
    dat_notland = dat_stack.where( np.isnan(landfrac_stack) | (landfrac_stack <= 0.5), drop=True)
    weights_notland = weights_stack.where( np.isnan(landfrac_stack) | (landfrac_stack <= 0.5), drop=True)

    dat_land_50n90n = dat_land.where( dat_land.lat < -50, drop=True)
    weights_land_50n90n = weights_land.where(weights_land.lat < -50, drop=True)

    dat_notland_50n90n = dat_notland.where( dat_notland.lat < -50, drop=True)
    weights_notland_50n90n = weights_notland.where(weights_notland.lat < -50, drop=True)

    dat_land_50n90n_w = dat_land_50n90n.weighted(weights_land_50n90n)
    dat_land_50n90n_m = dat_land_50n90n_w.sum("z")

    dat_notland_50n90n_w = dat_notland_50n90n.weighted(weights_notland_50n90n)
    dat_notland_50n90n_m = dat_notland_50n90n_w.sum("z")

    dat_land_50n90n_m = dat_land_50n90n_m.rename('land')
    dat_notland_50n90n_m = dat_notland_50n90n_m.rename('notland')

    dat = xr.merge([dat_land_50n90n_m, dat_notland_50n90n_m])

    dat.to_netcdf(pathout+'LENS2_'+varname+'_50s90s_land_notland.nc')


