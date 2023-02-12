import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial


from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig9/"

varnames=['AODVIS']

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds = ds.sel(year=slice(1920,2050))
    return ds

# 21 year running means
def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

dat =[]

for ivar in varnames:
    aer2 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
           "AAER_"+ivar+"_am.nc", preprocess=partial(preprocessor))

    aer2base = aer2.sel(year=slice(1920,1940)).mean('year')
    aer2anoms = aer2 - aer2base

    aer2gmt = avg.cosweightlonlat(aer2anoms[ivar], 0, 360, 50, 90)

    aer2gm = calc21ymean(aer2gmt)

    bootaer = boot.bootgen(aer2gm, nsamples = aer2gm.M.size, seed=4, nboots=1000)
    bootaer3 = boot.bootgen(aer2gm, nsamples = 3, seed=5, nboots=1000)

    bootaerm = bootaer.mean('isample')
    bootaer3m = bootaer3.mean('isample')

    min95 = bootaerm.quantile(0.025, dim='iboot') ; max95 = bootaerm.quantile(0.975, dim='iboot')
    min95_3 = bootaer3m.quantile(0.025, dim='iboot') ; max95_3 = bootaer3m.quantile(0.975, dim='iboot')

    min95 = min95.rename(ivar+'_min95') ; max95 = max95.rename(ivar+'_max95')
    min95_3 = min95_3.rename(ivar+'_min95_3') ; max95_3 = max95_3.rename(ivar+'_max95_3')

    if (ivar == varnames[0]):
        dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')
    else:
        dat = xr.merge([dat, min95, max95, min95_3, max95_3], compat='override')

dat.to_netcdf(pathout+'CESM2_AODVIS_bootstrap.nc')

