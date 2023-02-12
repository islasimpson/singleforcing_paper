import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

import dask
dask.config.set(**{'array.slicing.split_large_chunks': False})

from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig4/"

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds = ds.sel(year=slice(1920,2050))
    return ds

# 21 year running means
def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm


all_fsns = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
            "LENS1_FSNS_am.nc", preprocess=partial(preprocessor)).FSNS
xaer1_fsns = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
            "XAER_FSNS_am.nc", preprocess=partial(preprocessor)).FSNS

all_fsds = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
            "LENS1_FSDS_am.nc", preprocess=partial(preprocessor)).FSDS
xaer1_fsds = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
            "XAER_FSDS_am.nc", preprocess=partial(preprocessor)).FSDS

all_fsus = all_fsds - all_fsns
xaer1_fsus = xaer1_fsds - xaer1_fsns


all_fsns_gm = avg.cosweightlonlat(all_fsns, 0, 360, -90, 90)
xaer1_fsns_gm = avg.cosweightlonlat(xaer1_fsns, 0, 360, -90, 90)
all_fsds_gm = avg.cosweightlonlat(all_fsds, 0, 360, -90, 90)
xaer1_fsds_gm = avg.cosweightlonlat(xaer1_fsds, 0, 360, -90, 90)
all_fsus_gm = avg.cosweightlonlat(all_fsus, 0, 360, -90, 90)
xaer1_fsus_gm = avg.cosweightlonlat(xaer1_fsus, 0, 360, -90, 90)

all_alb = all_fsus_gm / all_fsds_gm
xaer1_alb = xaer1_fsus_gm / xaer1_fsds_gm

all_fsus_base = all_fsus_gm.sel(year=slice(1920,1940)).mean('year')
xaer1_fsus_base = xaer1_fsus_gm.sel(year=slice(1920,1940)).mean('year')

all_alb_base = all_alb.sel(year=slice(1920,1940)).mean('year')
xaer1_alb_base = xaer1_alb.sel(year=slice(1920,1940)).mean('year')

all_alb = all_alb - all_alb_base
xaer1_alb = xaer1_alb - xaer1_alb_base

all_fsus_gm = all_fsus_gm - all_fsus_base
xaer1_fsus_gm = xaer1_fsus_gm - xaer1_fsus_base

all_alb_21y = calc21ymean(all_alb)
xaer1_alb_21y = calc21ymean(xaer1_alb)
all_fsus_21y = calc21ymean(all_fsus_gm)
xaer1_fsus_21y = calc21ymean(xaer1_fsus_gm)

bootallalb = boot.bootgen(all_alb_21y, nsamples = all_alb_21y.M.size, seed=3, nboots=1000)
bootallfsus = boot.bootgen(all_fsus_21y, nsamples = all_fsus_21y.M.size, seed=4, nboots=1000)
bootxaeralb = boot.bootgen(xaer1_alb_21y, nsamples = xaer1_alb_21y.M.size, seed=3, nboots=1000)
bootxaerfsus = boot.bootgen(xaer1_fsus_21y, nsamples = xaer1_fsus_21y.M.size, seed=4, nboots=1000)
bootxaeralb3 = boot.bootgen(xaer1_alb_21y, nsamples = 3, seed=3, nboots=1000)
bootxaerfsus3 = boot.bootgen(xaer1_fsus_21y, nsamples = 3, seed=4, nboots=1000)

bootallalbm = bootallalb.mean('isample')
bootallfsusm = bootallfsus.mean('isample')
bootxaeralbm = bootxaeralb.mean('isample')
bootxaerfsusm = bootxaerfsus.mean('isample')
bootxaeralb3m = bootxaeralb3.mean('isample')
bootxaerfsus3m = bootxaerfsus3.mean('isample')

bootaeralb = bootallalbm - bootxaeralbm
bootaeralb3 = bootallalbm - bootxaeralb3m
bootaerfsus = bootallfsusm - bootxaerfsusm
bootaerfsus3 = bootallfsusm - bootxaerfsus3m

min95alb = bootaeralb.quantile(0.025, dim='iboot') ; max95alb = bootaeralb.quantile(0.975, dim='iboot')
min95alb3 = bootaeralb3.quantile(0.025, dim='iboot') ; max95alb3 = bootaeralb3.quantile(0.975, dim='iboot')
min95fsus = bootaerfsus.quantile(0.025, dim='iboot') ; max95fsus = bootaerfsus.quantile(0.975, dim='iboot')
min95fsus3 = bootaerfsus3.quantile(0.025, dim='iboot') ; max95fsus3 = bootaerfsus3.quantile(0.975, dim='iboot')

min95alb = min95alb.rename('Albedo_min95') ; max95alb = max95alb.rename('Albedo_max95')
min95alb3 = min95alb3.rename('Albedo_min95_3') ; max95alb3 = max95alb3.rename('Albedo_max95_3')
min95fsus = min95fsus.rename('FSUS_min95') ; max95fsus = max95fsus.rename('FSUS_max95')
min95fsus3 = min95fsus3.rename('FSUS_min95_3') ; max95fsus3 = max95fsus3.rename('FSUS_max95_3')


dat = xr.merge([min95alb, max95alb, min95alb3, max95alb3,
                min95fsus, max95fsus, min95fsus3, max95fsus3], compat='override')

dat.to_netcdf(pathout+'CESM1_bootstrap_Albedo_FSUS.nc')
