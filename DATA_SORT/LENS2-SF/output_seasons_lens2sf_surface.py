import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmemsGHG=15
nmemsAAER=15
nmemsBMB=15
nmemsEE=10 # !!! Update this when available

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/LENS2-SF/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX']
#varnames=['BURDENBCdn','BURDENDUSTdn','BURDENSO4dn','BURDENSOAdn']
#varnames=['TREFHT','AODVIS','BURDENBCdn','BURDENDUSTdn','BURDENPOMdn',
#          'BURDENSEASALTdn','BURDENSO4dn','BURDENSOAdn']

#varnames=['FSNS','FLNS','SHFLX','LHFLX','FLDS','FSDS',
#          'FLUT','FSNT','FSNSC','FSNTC']
varnames=['CLDTOT']


for varname in varnames:
    print(varname)
    topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/"+varname+"/"

    #----GHG
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsGHG+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-GHG."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-GHG-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+"GHG_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"GHG_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"GHG_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"GHG_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"GHG_"+varname+"_son.nc")
    #----End GHG

    #----AAER
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+"AAER_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"AAER_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"AAER_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"AAER_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"AAER_"+varname+"_son.nc")
    #----End AAER 

    #----BMB
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsBMB+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-BMB."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-BMB-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+"BMB_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"BMB_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"BMB_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"BMB_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"BMB_"+varname+"_son.nc")
    #----End BMB 

    #----EE
    memstr = [ str(100 + i) for i in np.arange(1,nmemsEE+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-EE."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-EE-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+"EE_"+varname+"_am.nc")
    djf.to_netcdf(pathout+"EE_"+varname+"_djf.nc")
    mam.to_netcdf(pathout+"EE_"+varname+"_mam.nc")
    jja.to_netcdf(pathout+"EE_"+varname+"_jja.nc")
    son.to_netcdf(pathout+"EE_"+varname+"_son.nc")
    #----EndEE
