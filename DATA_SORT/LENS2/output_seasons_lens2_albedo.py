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

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS2/"
topdir="/project/mojave/cesm2/LENS/atm/month_1/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX','FGR']
#varnames=['BURDENBCdn','BURDENDUSTdn','BURDENPOMdn','BURDENSEASALTdn','BURDENSO4dn','BURDENSOAdn' ]

memstr = lens.lens2memnamegen_second50(50)

varname='FSNS'
filelist=[sorted(glob.glob(topdir+varname+"/*.BHISTsmbb.*"+imem+"*.nc"))+
          sorted(glob.glob(topdir+varname+"/*.BSSP370smbb.*"+imem+"*.nc"))
          for imem in memstr ]
members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[[varname,'time_bnds']] for i in filelist]
fsns = xr.concat(members, dim='M', join='override', coords='minimal')
fsns = read.fixcesmtime(fsns)
fsns = fsns[varname]   

varname='FSDS'
filelist=[sorted(glob.glob(topdir+varname+"/*.BHISTsmbb.*"+imem+"*.nc"))+
          sorted(glob.glob(topdir+varname+"/*.BSSP370smbb.*"+imem+"*.nc"))
          for imem in memstr ]
members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[[varname,'time_bnds']] for i in filelist]
fsds = xr.concat(members, dim='M', join='override', coords='minimal')
fsds = read.fixcesmtime(fsds)
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
