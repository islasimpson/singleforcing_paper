import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

nmemsGHG=15
nmemsAAER=15
nmemsBMB=15
nmemsEE=15

topdir="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/RHO_top203/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"

#-------AAER
memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
filelist = [ sorted(glob.glob(topdir+"RHO_top203m_AAER2_"+imem+".nc")) for imem in memstr]
dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
dat = dat.where(dat.time.dt.month == 3, drop=True)
dat.load().to_netcdf(pathout+'AAER_March_RHO203.nc')




