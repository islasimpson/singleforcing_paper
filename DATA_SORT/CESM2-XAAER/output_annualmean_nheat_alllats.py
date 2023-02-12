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
topdir="/project/mojave/cesm2/Single_Forcing/ocn/month_1/N_HEAT/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"

memstr=[ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-xAER."+imem+".*")) for imem in memstr]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'], data_vars=['N_HEAT','time_bound'])
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time')
dat = dat.N_HEAT


heattransport45 = dat.isel(transport_reg=1).isel(transport_comp=[1,3,4]).sum('transport_comp')

heattransport45.load().to_netcdf(pathout+"xAER_N_HEAT_alllats_am.nc")

