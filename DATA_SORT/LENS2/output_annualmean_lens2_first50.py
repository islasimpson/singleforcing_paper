# Output annual means for the second 50 members of the CESM2 large ensemble
import xarray as xr
import numpy as np
import glob
import sys
import dask

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

import importlib
importlib.reload(lens)

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"
topdir="/project/mojave/cesm2/LENS/atm/month_1/"

#varnames=['AODVIS','FLNT','FSNT','FSNS','FSDS']
#varnames=['AODVIS','FSDS','FSNS']
varnames=['TREFHT']

memstr = lens.lens2memnamegen_first50(50)

for varname in varnames:
    print(varname)
    filelist=[sorted(glob.glob(topdir+varname+"/*.BHISTcmip6.*"+imem+"*.nc"))+
              sorted(glob.glob(topdir+varname+"/*.BSSP370cmip6.*"+imem+"*.nc"))
              for imem in memstr ]
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[[varname,'time_bnds']] for i in filelist]

    dat = xr.concat(members, dim='M', join='override', coords='minimal')
    dat = read.fixcesmtime(dat)
    dat = dat[varname]

    am = dat.groupby('time.year').mean('time').compute()
    am.to_netcdf(pathout+'LENS2_FIRST50_'+varname+'_am.nc')


