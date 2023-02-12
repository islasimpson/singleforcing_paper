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
#varnames=['TREFHT','FSNS','FSDS','ICEFRAC']
varnames=['FSNS','FSDS','TREFHT','ICEFRAC']

for varname in varnames:
    print(varname)
    topdir='/project/mojave/cesm2/Single_Forcing/atm/month_1/'+varname+'/'

    #------AAER
    memstr = [ str(i).zfill(3) for i in np.arange(1,nmemsAAER+1,1)]
    filelist = [ sorted(glob.glob(topdir+"*CESM2-SF-AAER."+imem+".*"))+
                 sorted(glob.glob(topdir+"*CESM2-SF-AAER-SSP370."+imem+"*")) for imem in memstr]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    dat = cal.season_ts(dat, 'DJF')
    dat.to_netcdf(pathout+'AAER_'+varname+'_djf.nc') 
