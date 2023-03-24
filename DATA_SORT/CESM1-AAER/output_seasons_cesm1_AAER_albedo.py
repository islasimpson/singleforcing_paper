import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/CESM1-AAER/"

#varnames=[ 'TREFHT','BURDENBC','BURDENDUST','BURDENPOM','BURDENSEASALT','BURDENSO4','BURDENSOA']
#varnames=[ 'FSNS','FLNS','SHFLX','LHFLX','FLDS','FSDS','FLUT','FSNT','FSNSC','FSNTC' ]

varname='FSDS'
topdir="/project/mojave/cesm1/Single_Forcing_asCESM2/atm/month_1/"+varname+"/"

memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
filelist = [ sorted(glob.glob(topdir+"*aaer*."+imem+"*.nc")) for imem in memstr ]
fsds = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
fsds = read.fixcesmtime(fsds)
fsds = fsds[varname]

varname='FSNS'
topdir="/project/mojave/cesm1/Single_Forcing_asCESM2/atm/month_1/"+varname+"/"

memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
filelist = [ sorted(glob.glob(topdir+"*aaer*."+imem+"*.nc")) for imem in memstr ]
fsns = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
fsns = read.fixcesmtime(fsns)
fsns = fsns[varname]

albedo = (fsds - fsns)/fsds

albedo = albedo.rename('Albedo')

am = albedo.groupby('time.year').mean('time').compute()
djf = cal.season_ts(albedo, 'DJF').compute()
mam = cal.season_ts(albedo, 'MAM').compute()
jja = cal.season_ts(albedo, 'JJA').compute()
son = cal.season_ts(albedo, 'SON').compute()

am.to_netcdf(pathout+'AAER_Albedo_am.nc')
djf.to_netcdf(pathout+'AAER_Albedo_djf.nc')
mam.to_netcdf(pathout+'AAER_Albedo_mam.nc')
jja.to_netcdf(pathout+'AAER_Albedo_jja.nc')
son.to_netcdf(pathout+'AAER_Albedo_son.nc')
