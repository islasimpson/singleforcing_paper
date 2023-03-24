# output seasonal means for the second 50 members of the CESM2 large ensemble
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

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS1/"
topdir="/project/mojave/cesm1/LENS/atm/month_1/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX','FGR']
#varnames=['AODVIS','BURDENBC','BURDENDUST','BURDENSEASALT','BURDENSO4','BURDENSOA']

memstr = lens.lens1memnamegen(40)

varname='FSNS'	
filelist = [sorted(glob.glob(topdir+varname+"/*.B20TRC5CNBDRD.f09_g16."+imem+"*.nc"))+
            sorted(glob.glob(topdir+varname+"/*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
            for imem in memstr]
members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], 
        coords='minimal')[[varname,'time_bnds']] for i in filelist]
members = [ read.fixcesmtime(members[i]).sel(time=slice("1920-01","2100-12")) 
         for i in np.arange(0,len(memstr),1)]

fsns = xr.concat(members, dim='M', join='override', coords='minimal')
fsns = fsns[varname]

varname='FSDS'	
filelist = [sorted(glob.glob(topdir+varname+"/*.B20TRC5CNBDRD.f09_g16."+imem+"*.nc"))+
            sorted(glob.glob(topdir+varname+"/*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
            for imem in memstr]
members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], 
        coords='minimal')[[varname,'time_bnds']] for i in filelist]
members = [ read.fixcesmtime(members[i]).sel(time=slice("1920-01","2100-12")) 
         for i in np.arange(0,len(memstr),1)]

fsds = xr.concat(members, dim='M', join='override', coords='minimal')
fsds = fsds[varname]

albedo = (fsds - fsns)/fsds

albedo = albedo.rename('Albedo')



am = albedo.groupby('time.year').mean('time').compute()
djf = cal.season_ts(albedo, "DJF").compute()
mam = cal.season_ts(albedo, "MAM").compute()
jja = cal.season_ts(albedo, "JJA").compute()
son = cal.season_ts(albedo, "SON").compute()

am.to_netcdf(pathout+"ALL_Albedo_am.nc")
djf.to_netcdf(pathout+"ALL_Albedo_djf.nc")
mam.to_netcdf(pathout+"ALL_Albedo_mam.nc")
jja.to_netcdf(pathout+"ALL_Albedo_jja.nc")
son.to_netcdf(pathout+"ALL_Albedo_son.nc")
