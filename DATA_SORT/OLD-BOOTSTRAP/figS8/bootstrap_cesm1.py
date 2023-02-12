import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

from math import nan

from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/figS8/"

varnames=['TREFHT','ICEFRAC']

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
notlandfrac = xr.where(np.isnan(landfrac.landfrac) | (landfrac.landfrac < 0.5), 1, nan)
landregions = xr.where(~np.isnan(landfrac.landfrac) | (landfrac.landfrac > 0.5), 1, nan)

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds['time'] = ds.time.dt.year
    ds = ds.rename({'time':'year'})
    ds = ds.sel(year=slice(1920,2050))
    return ds

def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

for ivar in varnames:

    xaer1 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
        "LENS1-SF/XAER_"+ivar+"_djf.nc", preprocess=partial(preprocessor))
    lens1 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
        "LENS1/LENS1_"+ivar+"_djf.nc", preprocess=partial(preprocessor))

    if (ivar == 'FSNO'):
        xaer1 = xaer1*landregions
        lens1 = lens1*landregions
    if (ivar == 'ICEFRAC'):
        xaer1 = xaer1*notlandfrac
        lens1 = lens1*notlandfrac

    xaer1base = xaer1.sel(year=slice(1920,1940)).mean('year')
    lens1base = lens1.sel(year=slice(1920,1940)).mean('year')

    xaer1 = xaer1 - xaer1base
    lens1 = lens1 - lens1base

    xaer1gmt = avg.cosweightlonlat(xaer1[ivar],0,360,-90,-50)
    xaer1gm = calc21ymean(xaer1gmt)

    lens1gmt = avg.cosweightlonlat(lens1[ivar],0,360,-90,-50)
    lens1gm = calc21ymean(lens1gmt)


    bootxaer1 = boot.bootgen(xaer1gm, nsamples = xaer1gm.M.size, seed=4, nboots=1000)
    bootxaer3 = boot.bootgen(xaer1gm, nsamples = 3, seed=5, nboots=1000)

    bootlens1 = boot.bootgen(lens1gm, nsamples = lens1gm.M.size, seed=5, nboots=1000)

    bootxaer1m = bootxaer1.mean('isample')
    bootxaer3m = bootxaer3.mean('isample')
    bootlens1m = bootlens1.mean('isample')

    bootaer1xway = bootlens1m - bootxaer1m
    bootaer1xway3 = bootlens1m - bootxaer3m

    min95 = bootaer1xway.quantile(0.025, dim='iboot') ; max95 = bootaer1xway.quantile(0.975, dim='iboot')
    min95_3 = bootaer1xway3.quantile(0.025, dim='iboot') ; max95_3 = bootaer1xway3.quantile(0.975, dim='iboot')

    min95 = min95.rename(ivar+'_min95') ; max95 = max95.rename(ivar+'_max95')
    min95_3 = min95_3.rename(ivar+'_min95_3') ; max95_3 = max95_3.rename(ivar+'_max95_3')

    if (ivar == varnames[0]):
        dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')
    else:
        dat = xr.merge([dat, min95, max95, min95_3, max95_3], compat='override')

dat.to_netcdf(pathout+'CESM1_bootstrap.nc')





