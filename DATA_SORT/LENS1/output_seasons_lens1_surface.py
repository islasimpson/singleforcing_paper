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

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS1/"
topdir="/project/mojave/cesm1/LENS/atm/month_1/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX','FGR']
#varnames=['AODVIS','BURDENBC','BURDENDUST','BURDENSEASALT','BURDENSO4','BURDENSOA']
varnames=['FLNTC']

memstr = lens.lens1memnamegen(40)

for varname in varnames:
    print(varname)
    filelist = [sorted(glob.glob(topdir+varname+"/*.B20TRC5CNBDRD.f09_g16."+imem+"*.nc"))+
                sorted(glob.glob(topdir+varname+"/*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
                for imem in memstr]
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'], 
            coords='minimal')[[varname,'time_bnds']] for i in filelist]
    members = [ read.fixcesmtime(members[i]).sel(time=slice("1920-01","2100-12")) 
             for i in np.arange(0,len(memstr),1)]

    dat = xr.concat(members, dim='M', join='override', coords='minimal')
    dat = dat[varname]

    am = dat.groupby('time.year').mean('time').compute()
#    djf = cal.season_ts(dat, "DJF").compute()
#    mam = cal.season_ts(dat, "MAM").compute()
#    jja = cal.season_ts(dat, "JJA").compute()
#    son = cal.season_ts(dat, "SON").compute()
#
    am.to_netcdf(pathout+"ALL_"+varname+"_am.nc")
#    djf.to_netcdf(pathout+"ALL_"+varname+"_djf.nc")
#    mam.to_netcdf(pathout+"ALL_"+varname+"_mam.nc")
#    jja.to_netcdf(pathout+"ALL_"+varname+"_jja.nc")
#    son.to_netcdf(pathout+"ALL_"+varname+"_son.nc")
