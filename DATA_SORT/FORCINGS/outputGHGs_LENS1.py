import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import readdata_utils as read
import glob

topdir="/project/mojave/cesm1/LENS/atm/month_1/PS/"
filelist =  sorted(glob.glob(topdir+'*B20TRC5CNBDRD.f09_g16.001*'))+\
             sorted(glob.glob(topdir+'*BRCP85C5CNBDRD.f09_g16.001*')) 
dat = xr.open_mfdataset(filelist)
dat = read.fixcesmtime(dat)
datam = dat.groupby('time.year').mean('time')

datam = datam[['co2vmr','ch4vmr','n2ovmr','f11vmr','f12vmr']]

datam.to_netcdf("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/GHGs_LENS1.nc")
