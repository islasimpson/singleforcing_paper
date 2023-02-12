import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/ocn/month_1/MOC/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/"

dat = xr.open_mfdataset(topdir+'*.MOC.1*.nc', decode_times='False')
timebnds = dat.time_bound
diff =  np.array(timebnds.isel(d2=1)) - np.array(timebnds.isel(d2=0))
diff = diff/2.
newtime = np.array(timebnds.isel(d2=0)) + diff
dat['time'] = newtime
dat = dat.groupby('time.year').mean('time').compute()
amoc = modes.calcAMOC(dat.MOC)
amoc.to_netcdf(pathout+'piControl_AMOC_am.nc')
