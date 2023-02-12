import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial


from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig3/"

varnames=['AODVIS','TREFHT','FSNT','FLNT','FSDS']

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
    all1 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
                "LENS1_"+ivar+"_am.nc", preprocess=partial(preprocessor))

    xaer1 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
                "XAER_"+ivar+"_am.nc", preprocess=partial(preprocessor))

    aer1 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"+\
               "AAER_"+ivar+"_am.nc", preprocess=partial(preprocessor))

    all1base = all1.sel(year=slice(1920,1940)).mean('year')
    xaer1base = xaer1.sel(year=slice(1920,1940)).mean('year')
    aer1base = aer1.sel(year=slice(1920,1940)).mean('year')

    all1 = all1 - all1base
    xaer1 = xaer1 - xaer1base
    aer1 = aer1 - aer1base

    all1gmt = avg.cosweightlonlat(all1[ivar], 0, 360, -90, 90)
    xaer1gmt = avg.cosweightlonlat(xaer1[ivar], 0, 360, -90, 90)
    aer1gmt = avg.cosweightlonlat(aer1[ivar], 0, 360, -90, 90)

    all1gm = calc21ymean(all1gmt)
    xaer1gm = calc21ymean(xaer1gmt)
    aer1gm = calc21ymean(aer1gmt)

    bootall = boot.bootgen(all1gm, nsamples = all1.M.size, seed=3, nboots=1000)
    bootxaer = boot.bootgen(xaer1gm, nsamples = xaer1gm.M.size, seed=4, nboots=1000)
    bootaer = boot.bootgen(aer1gm, nsamples = aer1gm.M.size, seed=5, nboots=1000)

    bootxaer3 = boot.bootgen(xaer1gm, nsamples = 3, seed=4, nboots=1000)

    bootallm = bootall.mean('isample')
    bootxaerm = bootxaer.mean('isample')
    bootxaer3m = bootxaer3.mean('isample')

    bootaerm = bootallm - bootxaerm
    bootaer3m = bootallm - bootxaer3m

    min95 = bootaerm.quantile(0.025, dim='iboot') ; max95 = bootaerm.quantile(0.975, dim='iboot')
    min95_3 = bootaer3m.quantile(0.025, dim='iboot') ; max95_3 = bootaer3m.quantile(0.975, dim='iboot')

    min95 = min95.rename(ivar+'_min95') ; max95 = max95.rename(ivar+'_max95')
    min95_3 = min95_3.rename(ivar+'_min95_3') ; max95_3 = max95_3.rename(ivar+'_max95_3')

    if (ivar == varnames[0]):
        dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')
    else:
        dat = xr.merge([dat, min95, max95, min95_3, max95_3], compat='override')


dat.to_netcdf(pathout+'CESM1_bootstrap.nc')
   

