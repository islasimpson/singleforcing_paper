import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

nmemsGHG=15
nmemsAAER=15
nmemsBMB=15
nmemsEE=15

topdir="/project/mojave/cesm2/Single_Forcing/ocn/month_1/MOC/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"

#------GHGs
#memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsGHG+1,1)]
#filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-GHG."+imem+".*"))+
#             sorted(glob.glob(topdir+"*CESM2-SF-GHG-SSP370."+imem+"*")) for imem in memstr]
#dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#dat = read.fixcesmtime(dat, timebndsvar='time_bound')
#dat = dat.groupby('time.year').mean('time').compute()
#amoc = modes.calcAMOC(dat.MOC)
#amoc.to_netcdf(pathout+'GHG_AMOC_am.nc')
#----------

#-------AAER
memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
             sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time').compute()
amoc = modes.calcAMOC(dat.MOC)
amoc.to_netcdf(pathout+'AAER_AMOC_am.nc')
#----------

#-------BMB
#memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsBMB+1,1)]
#filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-BMB."+imem+".*"))+
#             sorted(glob.glob(topdir+"*CESM2-SF-BMB-SSP370."+imem+"*")) for imem in memstr]
#dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#dat = read.fixcesmtime(dat, timebndsvar='time_bound')
#dat = dat.groupby('time.year').mean('time').compute()
#amoc = modes.calcAMOC(dat.MOC)
#amoc.to_netcdf(pathout+'BMB_AMOC_am.nc')
##----------
#
##-------EE
#memstr = [ str(i).zfill(3) for i in np.arange(101,101+nmemsEE+1,1)]
#filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-EE."+imem+".*"))+
#             sorted(glob.glob(topdir+"*CESM2-SF-EE-SSP370."+imem+"*")) for imem in memstr]
#dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#dat = read.fixcesmtime(dat, timebndsvar='time_bound')
#dat = dat.groupby('time.year').mean('time').compute()
#amoc = modes.calcAMOC(dat.MOC)
#amoc.to_netcdf(pathout+'EE_AMOC_am.nc')
