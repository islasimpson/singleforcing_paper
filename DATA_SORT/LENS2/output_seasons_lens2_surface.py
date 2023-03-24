# output seasonal means for the second 50 members of the CESM2 large ensemble
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

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS2/"
topdir="/project/mojave/cesm2/LENS/atm/month_1/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX','FGR']
#varnames=['BURDENBCdn','BURDENDUSTdn','BURDENPOMdn','BURDENSEASALTdn','BURDENSO4dn','BURDENSOAdn' ]

varnames=['SST']

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
 
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, "DJF").compute()
    mam = cal.season_ts(dat, "MAM").compute()
    jja = cal.season_ts(dat, "JJA").compute()
    son = cal.season_ts(dat, "SON").compute()

    am.to_netcdf(pathout+"ALL_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"ALL_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"ALL_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"ALL_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"ALL_"+varname+"_son.nc")
