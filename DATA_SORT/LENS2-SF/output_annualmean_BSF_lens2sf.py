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

topdir="/project/mojave/cesm2/Single_Forcing/ocn/month_1/BSF/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"

#-------AAER
memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
             sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'], data_vars=['BSF','time_bound'])
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.BSF
dat = dat.groupby('time.year').mean('time').compute()

#heattransport45 = dat.isel(transport_reg=1).isel(transport_comp=[1,3,4]).sum('transport_comp').sel(lat_aux_grid=50, method='nearest')


dat.load().to_netcdf(pathout+'AAER_BSF_am.nc')
#----------
