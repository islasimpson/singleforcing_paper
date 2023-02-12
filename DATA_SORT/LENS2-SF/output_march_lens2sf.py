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
nmemsEE=15

pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/'

#varnames=['AODVIS','FLNT','FSNT','FLNS','FSNS']
#varnames=['FSDS']
varnames=['PRECC','PRECL','SHFLX','QFLX']

for varname in varnames:
    print(varname)
    topdir='/project/mojave/cesm2/Single_Forcing/atm/month_1/'+varname+'/'

    #------GHGs
#    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsGHG+1,1)]
#    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-GHG."+imem+".*"))+
#                 sorted(glob.glob(topdir+"*CESM2-SF-GHG-SSP370."+imem+"*")) for imem in memstr]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'GHG_'+varname+'_am.nc')

    #------AAER
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    dat = dat.where(dat.time.dt.month == 3, drop=True)
    dat.to_netcdf(pathout+'AAER_'+varname+'_march.nc')
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'AAER_'+varname+'_am.nc') 

    #------BMB
#    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsBMB+1,1)]
#    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-BMB."+imem+".*"))+
#                 sorted(glob.glob(topdir+"*CESM2-SF-BMB-SSP370."+imem+"*")) for imem in memstr]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'BMB_'+varname+'_am.nc')

    #------EE
#    memstr = [ str(i).zfill(3) for i in np.arange(101,101+nmemsEE+1,1)]
#    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-EE."+imem+".*"))+
#                 sorted(glob.glob(topdir+"*CESM2-SF-EE-SSP370."+imem+"*")) for imem in memstr]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'EE_'+varname+'_am.nc')


