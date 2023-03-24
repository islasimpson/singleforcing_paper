import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing/DATA_SORT/CESM2-XAER/"

#varnames=['TREFHT','FLNS','FSNS','SHFLX','LHFLX']
#varnames=['BURDENBCdn','BURDENDUSTdn','BURDENPOMdn','BURDENSO4dn','BURDENSOAdn']
#varnames=['AODdnDUST1','AODdnDUST2','AODdnDUST3','AODdnMODE1','AODdnMODE2','AODdnMODE3',
#    'AODDUST','AODDUST1','AODDUST2','AODDUST3','AODPOMdn','AODSO4dn','AODSOAdn','AODSSdn']
varnames=['CLDTOT']

for varname in varnames:
    print(varname)
    topdir="/project/mojave/cesm2/Single_Forcing/atm/month_1/"+varname+"/"

    memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-xAER."+imem+".*")) for imem in memstr ]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    am = dat.groupby('time.year').mean('time').compute()
    djf = cal.season_ts(dat, 'DJF').compute()
    mam = cal.season_ts(dat, 'MAM').compute()
    jja = cal.season_ts(dat, 'JJA').compute()
    son = cal.season_ts(dat, 'SON').compute()

    am.to_netcdf(pathout+'xAER_'+varname+'_am.nc')
    djf.to_netcdf(pathout+'xAER_'+varname+'_djf.nc')
    mam.to_netcdf(pathout+'xAER_'+varname+'_mam.nc')
    jja.to_netcdf(pathout+'xAER_'+varname+'_jja.nc')
    son.to_netcdf(pathout+'xAER_'+varname+'_son.nc')
