import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm1/LENS/ocn/month_1/HMXL/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"

memstr = lens.lens1memnamegen(40)

filelist = [ sorted(glob.glob(topdir+"*.B20TRC5CNBDRD.f09_g16."+imem+"*.nc"))+
             sorted(glob.glob(topdir+"*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
             for imem in memstr ]


members = [ xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal', compat='override')[['HMXL','time_bound']].sel(time=slice("1920-02","2100-01")) for i in filelist ]

#members = [ xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal', compat='override')[['HMXL','time_bound']].sel(time=slice("1920-02","2100-01"))
#            for i in filelist ]

dat = xr.concat(members, dim='M', join='override', coords='minimal', compat='override')
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.where(dat.time.dt.month == 3, drop=True)
dat.load().to_netcdf(pathout+'LENS1_March_HMXL.nc')
