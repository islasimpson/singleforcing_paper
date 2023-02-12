import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from CASutils import readdata_utils as read
import glob

topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/SOLIN/"
memstr = [ str(imem).zfill(3) for imem in np.arange(101,115+1,1) ]
filelist = [ sorted(glob.glob(topdir+'*CESM2-SF-EE.'+imem+'*'))+\
             sorted(glob.glob(topdir+'*CESM2-SF-EE-SSP370.'+imem+'*')) for imem in memstr ]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
dat = read.fixcesmtime(dat)
datam = dat.groupby('time.year').mean('time')

datam = datam.SOLIN

datam.to_netcdf("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FORCINGS/SOLIN.nc")
