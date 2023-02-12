import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/atm/month_1/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/"

fsns = xr.open_mfdataset(topdir+"FSNS/*.nc", decode_times='False', data_vars=['FSNS','time_bnds'])
timebnds = fsns.time_bnds
diff =  np.array(timebnds.isel(nbnd=1)) - np.array(timebnds.isel(nbnd=0))
diff = diff/2.
newtime = np.array(timebnds.isel(nbnd=0)) + diff
fsns['time'] = newtime
fsnsam = fsns.groupby('time.year').mean('time')


fsds = xr.open_mfdataset(topdir+"FSDS/*.nc", decode_times='False', data_vars=['FSDS','time_bnds'])
timebnds = fsds.time_bnds
diff =  np.array(timebnds.isel(nbnd=1)) - np.array(timebnds.isel(nbnd=0))
diff = diff/2.
newtime = np.array(timebnds.isel(nbnd=0)) + diff
fsds['time'] = newtime
fsdsam = fsds.groupby('time.year').mean('time')

fsusam = fsdsam.FSDS - fsnsam.FSNS
del(fsnsam)
fsusam = fsusam.rename('FSUS')

fsdsgm = avg.cosweightlonlat(fsdsam, 0, 360, -90, 90)
fsusgm = avg.cosweightlonlat(fsusam, 0, 360, -90, 90)

albedo = fsusgm / fsdsgm.FSDS

albedo = albedo.rename('Albedo')

albedo.to_netcdf(pathout+"piControl_albedo_am.nc")
