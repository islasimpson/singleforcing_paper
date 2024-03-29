import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/CESM1-AAER/"

#varnames=[ 'TREFHT','BURDENBC','BURDENDUST','BURDENPOM','BURDENSEASALT','BURDENSO4','BURDENSOA']
#varnames=[ 'FSNS','FLNS','SHFLX','LHFLX','FLDS','FSDS','FLUT','FSNT','FSNSC','FSNTC' ]
varnames=['FLNTC']

for varname in varnames:
    print(varname)
    topdir="/project/mojave/cesm1/Single_Forcing_asCESM2/atm/month_1/"+varname+"/"

    memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*aaer*."+imem+"*.nc")) for imem in memstr ]

    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+'AAER_'+varname+'_am.nc')
    djf.to_netcdf(pathout+'AAER_'+varname+'_djf.nc')
    mam.to_netcdf(pathout+'AAER_'+varname+'_mam.nc')
    jja.to_netcdf(pathout+'AAER_'+varname+'_jja.nc')
    son.to_netcdf(pathout+'AAER_'+varname+'_son.nc')
