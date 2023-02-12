import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'
topdir='/project/mojave/cesm1/Single_Forcing/atm/month_1/'

nmemsxghg=20
nmemsxaer=20
nmemsxbmb=15

#varnames=['AODVIS','FLNT','FSNT','FLNS','FSNS']
#varnames=['FSNTC','FLNTC']
varnames=['FSNO']

for varname in varnames:
    print(varname)
    topdir='/project/mojave/cesm1/Single_Forcing/lnd/month_1/'+varname+'/'

    #-----XGHG
#    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxghg+1,1)]
#    filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xghg.'+imem+'.*'))
#                 for imem in memstr ]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'XGHG_'+varname+'_am.nc')

    #-----XAER
    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
    filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
                 for imem in memstr ]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'], data_vars=[varname,'time_bounds'], compat='override', coords='minimal')
    dat = read.fixcesmtime(dat, timebndsvar='time_bounds')
    dat = dat[varname]

    jja = cal.season_ts(dat, 'JJA')

    jja.to_netcdf(pathout+'XAER_'+varname+'_jja.nc')

    #-----XBMB
#    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxbmb+1,1)]
#    filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xbmb.'+imem+'.*'))
#                 for imem in memstr ]
#    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
#    dat = read.fixcesmtime(dat)
#    dat = dat[varname]
#    am = dat.groupby('time.year').mean('time').compute()
#    am.to_netcdf(pathout+'XBMB_'+varname+'_am.nc')
