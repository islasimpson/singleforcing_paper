import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import modes_utils as modes

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"

topdir="/project/mojave/cesm1/Single_Forcing_asCESM2/ocn/month_1/MOC/"

memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
filelist = [ sorted(glob.glob(topdir+'*aaer*.'+imem+'*.nc')) for imem in memstr]

dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
#dat = dat[varname]
dat = dat.groupby('time.year').mean('time').compute()
amoc = modes.calcAMOC(dat.MOC, lat=45.)

amoc.load().to_netcdf(pathout+'AAER_amoc_45_am.nc')
