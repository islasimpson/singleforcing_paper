import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"

#varnames=['AODVIS','TREFHT','FSNT','FLNT','FSNS','FLNS']
varnames=['FSNO']

for varname in varnames:
    print(varname)
    topdir="/project/mojave/cesm2/Single_Forcing/lnd/month_1/"+varname+"/"

    memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-xAER."+imem+".*")) for imem in memstr ]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'],
           data_vars=[varname,'time_bounds'], compat='override', coords='minimal')
    dat = read.fixcesmtime(dat, timebndsvar='time_bounds')
    dat = dat[varname]
    jja = cal.season_ts(dat,'JJA')

    jja.to_netcdf(pathout+'xAER_'+varname+'_jja.nc')
