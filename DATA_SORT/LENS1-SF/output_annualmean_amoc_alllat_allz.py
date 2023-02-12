import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import modes_utils as modes

pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'
topdir='/project/mojave/cesm1/Single_Forcing/atm/month_1/'

nmemsxghg=20
nmemsxaer=20
nmemsxbmb=15

topdir="/project/mojave/cesm1/Single_Forcing/ocn/month_1/MOC/"
pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'

memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
             for imem in memstr ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])[['MOC','time_bound']]
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time')
amoc = dat.MOC.sum('moc_comp').isel(transport_reg=1)

amoc.load().to_netcdf(pathout+'XAER_AMOC_am_alllat_allz.nc')
