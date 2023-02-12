import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

nmemsxaer=20

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"

vars=( "SFWF", "QFLUX", "SHF" )

for ivar in vars:
    print(ivar)
    topdir="/project/mojave/cesm1/Single_Forcing/ocn/month_1/"+ivar+"/"
    memstr = [str(i).zfill(3) for i in np.arange(1,nmemsxaer+1,1)]
    filelist = [ sorted(glob.glob(topdir+'b.e11.B20TRLENS_RCP85.f09_g16.xaer.'+imem+'.*'))
             for imem in memstr ]
    dat = xr.open_mfdataset(filelist, combine='nested', concat_dim=['M','time'])[[ivar,'time_bound']]
    dat = read.fixcesmtime(dat, timebndsvar='time_bound')
    dat = dat[ivar]
    dat = cal.season_ts(dat, season='DJFM')
    dat.to_netcdf(pathout+'xAER_'+ivar+'_DJFM.nc')
