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
varnames=['FSDS','FSNS']

for varname in varnames:
    print(varname)
    topdir='/project/mojave/cesm1/Single_Forcing/atm/month_1/'+varname+'/'

    #-----XAER
    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
    filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
                 for imem in memstr ]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])
    dat = read.fixcesmtime(dat)
    dat = dat[varname]
    dat= dat.mean('M').compute()
    dat.to_netcdf(pathout+'XAER_'+varname+'_ensemblemean_allmonths.nc')
