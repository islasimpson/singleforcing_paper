import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial


from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig3/"

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds = ds.sel(year=slice(1920,2050))
    return ds

# 21 year running means
def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm


all_fsnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
            "LENS1_FSNT_am.nc", preprocess=partial(preprocessor)).FSNT
xaer1_fsnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
            "XAER_FSNT_am.nc", preprocess=partial(preprocessor)).FSNT

all_flnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
            "LENS1_FLNT_am.nc", preprocess=partial(preprocessor)).FLNT
xaer1_flnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
            "XAER_FLNT_am.nc", preprocess=partial(preprocessor)).FLNT

all_fsnt_base = all_fsnt.sel(year=slice(1920,1940)).mean('year')
xaer1_fsnt_base = xaer1_fsnt.sel(year=slice(1920,1940)).mean('year')

all_flnt_base = all_flnt.sel(year=slice(1920,1940)).mean('year')
xaer1_flnt_base = xaer1_flnt.sel(year=slice(1920,1940)).mean('year')

all_fsnt = all_fsnt - all_fsnt_base
all_flnt = all_flnt - all_flnt_base
xaer1_fsnt = xaer1_fsnt - xaer1_fsnt_base
xaer1_flnt = xaer1_flnt - xaer1_flnt_base

all_fsnt_gmt = avg.cosweightlonlat(all_fsnt, 0, 360, -90, 90)
all_flnt_gmt = avg.cosweightlonlat(all_flnt, 0, 360, -90, 90)
xaer1_fsnt_gmt = avg.cosweightlonlat(xaer1_fsnt, 0, 360, -90, 90)
xaer1_flnt_gmt = avg.cosweightlonlat(xaer1_flnt, 0, 360, -90, 90)

all_fsnt_gm = calc21ymean(all_fsnt_gmt)
all_flnt_gm = calc21ymean(all_flnt_gmt)
xaer1_fsnt_gm = calc21ymean(xaer1_fsnt_gmt)
xaer1_flnt_gm = calc21ymean(xaer1_flnt_gmt)

all_toa = all_fsnt_gm - all_flnt_gm
xaer1_toa = xaer1_fsnt_gm - xaer1_flnt_gm

bootall = boot.bootgen(all_toa, nsamples = all_toa.M.size, seed=3, nboots=1000)
bootxaer = boot.bootgen(xaer1_toa, nsamples = xaer1_toa.M.size, seed=4, nboots=1000)
bootxaer3 = boot.bootgen(xaer1_toa, nsamples = 3, seed=5, nboots=1000)

bootallm = bootall.mean('isample')
bootxaerm = bootxaer.mean('isample')
bootxaer3m = bootxaer3.mean('isample')

bootaerm = bootallm - bootxaerm
bootaer3m = bootallm - bootxaer3m

min95 = bootaerm.quantile(0.025, dim='iboot') ; max95 = bootaerm.quantile(0.975, dim='iboot')
min95_3 = bootaer3m.quantile(0.025, dim='iboot') ; max95_3 = bootaer3m.quantile(0.975, dim='iboot')

min95 = min95.rename('TOA_min95') ; max95 = max95.rename('TOA_max95')
min95_3 = min95_3.rename('TOA_min95_3') ; max95_3 = max95_3.rename('TOA_max95_3')

dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')

dat.to_netcdf(pathout+'CESM1_bootstrap_TOA.nc')
