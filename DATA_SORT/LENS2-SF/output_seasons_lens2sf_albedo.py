import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmemsGHG=15
nmemsAAER=15
nmemsBMB=15
nmemsEE=10 # !!! Update this when available

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS2-SF/"


#----AAER
varname='FSDS'
topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/"+varname+"/"
memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
             sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
fsds = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
fsds = read.fixcesmtime(fsds)
fsds = fsds[varname]

varname='FSNS'
topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/"+varname+"/"
memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
             sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
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

am.to_netcdf(pathout+"AAER_Albedo_am.nc")
djf.to_netcdf(pathout+"AAER_Albedo_djf.nc")
mam.to_netcdf(pathout+"AAER_Albedo_mam.nc")
jja.to_netcdf(pathout+"AAER_Albedo_jja.nc")
son.to_netcdf(pathout+"AAER_Albedo_son.nc")
#----End AAER 
