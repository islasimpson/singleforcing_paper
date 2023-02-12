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
nmemsEE=15

pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/'

#varnames=['AODVIS','FLNT','FSNT','FLNS','FSNS']
#varnames=['FSDS']
varnames=['FSNO']

for varname in varnames:
    print(varname)
    topdir='/project/mojave/cesm2/Single_Forcing/lnd/month_1/'+varname+'/'

    #------AAER
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'],
             data_vars=[varname,'time_bounds'], compat='override', coords='minimal')
    dat = read.fixcesmtime(dat, timebndsvar='time_bounds')
    dat = dat[varname]
    dat = cal.season_ts(dat, 'JJA')
    dat.to_netcdf(pathout+'AAER_'+varname+'_jja.nc') 
