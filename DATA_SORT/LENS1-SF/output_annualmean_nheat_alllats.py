import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal
from CASutils import modes_utils as modes

nmemsxghg=20
nmemsxaer=20
nmemsxbmb=15

topdir="/project/mojave/cesm1/Single_Forcing/ocn/month_1/N_HEAT/"
pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'

memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
             for imem in memstr ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])[['N_HEAT','time_bound']]
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time')

heattransport = dat.isel(transport_reg=1).isel(transport_comp=[1,3,4]).sum('transport_comp')

heattransport.load().to_netcdf(pathout+'XAER_N_HEAT_alllats_am.nc')
