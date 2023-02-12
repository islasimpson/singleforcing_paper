import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

topdir="/project/mojave/cesm1/LENS/atm/month_1/CLDTOT/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM1/"

trefht = xr.open_mfdataset(topdir+'b.e11.B1850C5CN.f09_g16.005*.nc', 
           decode_times='False', data_vars=['CLDTOT', 'time_bnds'])
timebnds = trefht.time_bnds
diff =  np.array(timebnds.isel(nbnd=1)) - np.array(timebnds.isel(nbnd=0))
diff = diff/2.
newtime = np.array(timebnds.isel(nbnd=0)) + diff
trefht['time'] = newtime
trefhtam = trefht.groupby('time.year').mean('time').load()

trefhtam.to_netcdf(pathout+'picontrol_CLDTOT_am.nc')
