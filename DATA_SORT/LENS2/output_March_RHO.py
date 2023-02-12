import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

import dask
dask.config.set(**{'array.slicing.split_large_chunks': True})

topdir="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/RHO_top203m/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"

memstr = lens.lens2memnamegen_second50(50)

filelist = [ sorted(glob.glob(topdir+"*"+imem+"*.nc")) for imem in memstr ]
members = [ xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal') for i in filelist ]

dat = xr.concat(members, dim='M', join='override', coords='minimal')
dat = dat.where(dat.time.dt.month == 3, drop=True)
dat.load().to_netcdf(pathout+"LENS2_March_RHO203.nc")
