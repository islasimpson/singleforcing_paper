import xarray as xr
import numpy as np
import glob
import sys
import pandas as pd

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import readdata_utils as read

topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/atm/month_1/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/"

varname=['TREFHT']

for ivar in varname:
    path = topdir+ivar+'/'
    dat = xr.open_mfdataset(topdir+ivar+"/*."+ivar+".1*.nc", coords='minimal', decode_times='False')
    timebnds = dat.time_bnds
    diff = np.array(timebnds.isel(nbnd=1)) - np.array(timebnds.isel(nbnd=0))
    diff = diff/2.
    newtime = np.array(timebnds.isel(nbnd=0)) + diff
    dat['time'] = newtime
    dat = dat.groupby('time.year').mean('time')
    dat.to_netcdf(pathout+"TREFHT_picontrol_am.nc")
