import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm1/Single_Forcing/ocn/month_1/HMXL/"
pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'

nmemsxaer=20

memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
             for imem in memstr ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])[['HMXL','time_bound']]
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.where(dat.time.dt.month == 3, drop=True)
dat.load().to_netcdf(pathout+'XAER_March_HMXL.nc')

