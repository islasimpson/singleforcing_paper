import xarray as xr
import numpy as np
import glob
import sys

from CASutils import averaging_utils as avg
from CASutils import readdata_utils as read
from CASutils import calendar_utils as cal

nmems=3

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"

#varnames=['AODVIS','TREFHT','FLNT','FSNT','FLNS','FSNS']
#varnames=['TREFHT','FSNS','FSDS','ICEFRAC']
varnames=['FSNO']


for varname in varnames:
    print(varname)
    topdir="/project/mojave/cesm1/Single_Forcing_asCESM2/lnd/month_1/"+varname+"/"

    memstr = [ str(i).zfill(3) for i in np.arange(1,nmems+1,1)]
    filelist = [ sorted(glob.glob(topdir+'*aaer*.'+imem+'*.nc')) for imem in memstr]

    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat, timebndsvar='time_bounds')
    dat = dat[varname]

    jja = cal.season_ts(dat, 'JJA')

    jja.to_netcdf(pathout+'AAER_'+varname+'_jja.nc')
