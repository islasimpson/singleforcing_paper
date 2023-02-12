import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/atm/month_1/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/"

trefht = xr.open_mfdataset(topdir+"BURDENSO4dn/*190001-200012.nc", decode_times='False', data_vars=['BURDENSO4dn','time_bnds'])
timebnds = trefht.time_bnds
diff =  np.array(timebnds.isel(nbnd=1)) - np.array(timebnds.isel(nbnd=0))
diff = diff/2.
newtime = np.array(timebnds.isel(nbnd=0)) + diff
trefht['time'] = newtime
trefhtam = trefht.groupby('time.year').mean('time')

trefhtam.to_netcdf(pathout+'picontrol_BURDENSO4dn_am.nc')

