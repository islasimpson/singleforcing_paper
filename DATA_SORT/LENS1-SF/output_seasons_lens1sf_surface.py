import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS1-SF/"
topdir="/project/mojave/cesm1/Single_Forcing/atm/month_1/"

nmemsxghg=20
nmemsxaer=20
nmemsxbmb=15

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX']
#varnames=['BURDENBC','BURDENDUST','BURDENSEASALT','BURDENSO4','BURDENSOA']
#varnames=['FLUT','FSNT','FLDS']
varnames=['CLDTOT']

for varname in varnames:
    print(varname) 
    topdir="/project/mojave/cesm1/Single_Forcing/atm/month_1/"+varname+"/"

    #----------------XGHG
#    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxghg+1,1)]
#    filelist = [sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xghg.'+imem+'.*'))
#                   for imem in memstr]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    djf = cal.season_ts(dat, 'DJF').compute()
#    mam = cal.season_ts(dat, 'MAM').compute()
#    jja = cal.season_ts(dat, 'JJA').compute()
#    son = cal.season_ts(dat, 'SON').compute()

#    am.to_netcdf(pathout+"XGHG_"+varname+"_am.nc")
#    djf.to_netcdf(pathout+"XGHG_"+varname+"_djf.nc")
#    mam.to_netcdf(pathout+"XGHG_"+varname+"_mam.nc")
#    jja.to_netcdf(pathout+"XGHG_"+varname+"_jja.nc")
#    son.to_netcdf(pathout+"XGHG_"+varname+"_son.nc")
    #-------End XGHG

    #----------------XAER
    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
    filelist = [sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
                   for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
#    djf = cal.season_ts(dat, 'DJF').compute()
#    mam = cal.season_ts(dat, 'MAM').compute()
#    jja = cal.season_ts(dat, 'JJA').compute()
#    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+"XAER_"+varname+"_am.nc")
#    djf.to_netcdf(pathout+"XAER_"+varname+"_djf.nc")
#    mam.to_netcdf(pathout+"XAER_"+varname+"_mam.nc")
#    jja.to_netcdf(pathout+"XAER_"+varname+"_jja.nc")
#    son.to_netcdf(pathout+"XAER_"+varname+"_son.nc")
    #-------End XGHG

    #----------------XBMB
#    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxbmb+1,1)]
#    filelist = [sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xbmb.'+imem+'.*'))
#                   for imem in memstr]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    djf = cal.season_ts(dat, 'DJF').compute()
#    mam = cal.season_ts(dat, 'MAM').compute()
#    jja = cal.season_ts(dat, 'JJA').compute()
#    son = cal.season_ts(dat, 'SON').compute()

#    am.to_netcdf(pathout+"XBMB_"+varname+"_am.nc")
#    djf.to_netcdf(pathout+"XBMB_"+varname+"_djf.nc")
#    mam.to_netcdf(pathout+"XBMB_"+varname+"_mam.nc")
#    jja.to_netcdf(pathout+"XBMB_"+varname+"_jja.nc")
#    son.to_netcdf(pathout+"XBMB_"+varname+"_son.nc")
    #-------End XGHG




