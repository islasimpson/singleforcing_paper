import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS1-SF/"
topdir="/project/mojave/cesm1/Single_Forcing/atm/month_1/"

nmemsxghg=20
nmemsxaer=20
nmemsxbmb=15

#----------------XAER
varname='FSNS'
topdir="/project/mojave/cesm1/Single_Forcing/atm/month_1/"+varname+"/"
memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
filelist = [sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
               for imem in memstr]
fsns = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
fsns = read.fixcesmtime(fsns)
fsns = fsns[varname]

varname='FSDS'
topdir="/project/mojave/cesm1/Single_Forcing/atm/month_1/"+varname+"/"
memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
filelist = [sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
               for imem in memstr]
fsds = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
fsds = read.fixcesmtime(fsds)
fsds = fsds[varname]

albedo = (fsds - fsns)/fsds
albedo = albedo.rename('Albedo')

am = albedo.groupby('time.year').mean('time').compute()
djf = cal.season_ts(albedo, 'DJF').compute()
mam = cal.season_ts(albedo, 'MAM').compute()
jja = cal.season_ts(albedo, 'JJA').compute()
son = cal.season_ts(albedo, 'SON').compute()

am.to_netcdf(pathout+"XAER_Albedo_am.nc")
djf.to_netcdf(pathout+"XAER_Albedo_djf.nc")
mam.to_netcdf(pathout+"XAER_Albedo_mam.nc")
jja.to_netcdf(pathout+"XAER_Albedo_jja.nc")
son.to_netcdf(pathout+"XAER_Albedo_son.nc")
