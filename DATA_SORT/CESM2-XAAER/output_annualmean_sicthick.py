import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

nmems=3
topdir="/project/mojave/cesm2/Single_Forcing/ice/month_1/hi/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"

memstr=[ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-xAER."+imem+".*")) for imem in memstr]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'], data_vars=['hi','time_bounds'])
dat = read.fixcesmtime(dat, timebndsvar='time_bounds')
dat = dat.hi
dat = dat.groupby('time.year').mean('time')
dat.load().to_netcdf(pathout+'xAER_sicthick_am.nc')
