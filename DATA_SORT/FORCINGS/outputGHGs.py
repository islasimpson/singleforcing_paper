import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import readdata_utils as read
import glob

topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/PS/"
filelist =  sorted(glob.glob(topdir+'*CESM2-SF-GHG.001*'))+\
             sorted(glob.glob(topdir+'*CESM2-SF-GHG-SSP370.001*')) 
dat = xr.open_mfdataset(filelist)
dat = read.fixcesmtime(dat)
datam = dat.groupby('time.year').mean('time')

datam = datam[['co2vmr','ch4vmr','n2ovmr','f11vmr','f12vmr']]

datam.to_netcdf("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/GHGs.nc")
