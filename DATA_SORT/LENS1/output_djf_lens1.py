# output annual means for the CESM1 large ensemble
import xarray as xr
import numpy as np
import glob
import sys
import dask

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

pathout='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/'
topdir='/project/mojave/cesm1/LENS/atm/month_1/'

#varnames=['AODVIS','FLNT','FSNT','FLNS','FSNS']
#varnames=['FSDS']
varnames=['TREFHT','FSNS','FSDS','ICEFRAC']

memstr = lens.lens1memnamegen(40)

for varname in varnames:
    print(varname)
    # Get the filenames
    filelist = [sorted(glob.glob(topdir+varname+"/*.B20TRC5CNBDRD.f09_g16."+imem+"*.nc"))+
                sorted(glob.glob(topdir+varname+"/*.BRCP85C5CNBDRD.f09_g16."+imem+"*.nc"))
                for imem in memstr]
    # read in each member
    members = [xr.open_mfdataset(i, combine='nested', concat_dim=['time'],
            coords='minimal')[[varname,'time_bnds']] for i in filelist]
    # Fix the CESM time axis
    members = [ read.fixcesmtime(members[i]).sel(time=slice("1920-01","2100-12"))
             for i in np.arange(0,len(memstr),1)]

    # Concatenate along the member axis
    dat = xr.concat(members, dim='M', join='override', coords='minimal')
    dat = dat[varname]

    jja = cal.season_ts(dat, 'DJF')

    jja.to_netcdf(pathout+'LENS1_'+varname+'_djf.nc')
