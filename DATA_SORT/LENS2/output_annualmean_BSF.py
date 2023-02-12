import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm2/LENS/ocn/month_1/BSF/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"

memstr = lens.lens2memnamegen_second50(50)

filelist = [ sorted(glob.glob(topdir+"*.BHISTsmbb.*"+imem+"*.nc"))+
             sorted(glob.glob(topdir+"*.BSSP370smbb.*"+imem+"*.nc"))
             for imem in memstr ] 
members = [ xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[['BSF','time_bound']]
            for i in filelist ]
dat = xr.concat(members, dim='M', join='override', coords='minimal')
dat = read.fixcesmtime(dat, timebndsvar='time_bound')
dat = dat.groupby('time.year').mean('time')

#heattransport45 = dat.isel(transport_reg=1).isel(transport_comp=[1,3,4]).sum('transport_comp').sel(lat_aux_grid=50, method='nearest')

dat.to_netcdf(pathout+'LENS2_BSF_am.nc')
