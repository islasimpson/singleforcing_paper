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
#varnames=['TREFHT','FSNS','FSDS','ICEFRAC']
varnames=['FSNS','FSDS']

memstr = lens.lens2memnamegen_second50(50)

for varname in varnames:
    print(varname)
    filelist=[sorted(glob.glob(topdir+varname+"/*.BHISTsmbb.*"+imem+"*.nc"))+
              sorted(glob.glob(topdir+varname+"/*.BSSP370smbb.*"+imem+"*.nc"))
              for imem in memstr ]
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], coords='minimal')[[varname,'time_bnds']] for i in filelist]

    dat = xr.concat(members, dim='M', join='override', coords='minimal')
    dat = read.fixcesmtime(dat)
    dat = dat[varname]

    jja = cal.season_ts(dat, 'JJA')

    jja.to_netcdf(pathout+'LENS2_'+varname+'_jja.nc')


