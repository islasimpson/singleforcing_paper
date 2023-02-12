import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm1/LENS/ocn/month_1/MOC/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"

memstr = lens.lens1memnamegen(40)

filelist = [ sorted(glob.glob(topdir+"*.B5CNBDRD.f09_g16."+imem+"*.nc"))+
             sorted(glob.glob(topdir+"*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
             for imem in memstr ]
members = [ xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[['MOC','time_bound']]
            for i in filelist ]
dat = xr.concat(members, dim='M', join='override', coords='minimal')
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time')
amoc = dat.MOC.sum('moc_comp').isel(transport_reg=1)
#amoc = modes.calcAMOC(dat.MOC, lat=45.)

amoc.load().to_netcdf(pathout+'LENS1_AMOC_am_alllat_allz.nc')

