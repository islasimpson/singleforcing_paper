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

aer2_fsnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
           "AAER_FSNT_am.nc", preprocess=partial(preprocessor)).FSNT
aer2_flnt = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
           "AAER_FLNT_am.nc", preprocess=partial(preprocessor)).FLNT

aer2_fsnt_base = aer2_fsnt.sel(year=slice(1920,1940)).mean('year')
aer2_flnt_base = aer2_flnt.sel(year=slice(1920,1940)).mean('year')

aer2_fsnt = aer2_fsnt - aer2_fsnt_base
aer2_flnt = aer2_flnt - aer2_flnt_base

aer2_fsnt_gmt = avg.cosweightlonlat(aer2_fsnt, 0, 360, -90, 90)
aer2_flnt_gmt = avg.cosweightlonlat(aer2_flnt, 0, 360, -90, 90)

aer2_fsnt_gm = calc21ymean(aer2_fsnt_gmt)
aer2_flnt_gm = calc21ymean(aer2_flnt_gmt)

aer2_toa = aer2_fsnt_gm - aer2_flnt_gm
    
bootaer = boot.bootgen(aer2_toa, nsamples = aer2_toa.M.size, seed=4, nboots=1000)
bootaer3 = boot.bootgen(aer2_toa, nsamples = 3, seed=5, nboots=1000)

bootaerm = bootaer.mean('isample')
bootaer3m = bootaer3.mean('isample')

min95 = bootaerm.quantile(0.025, dim='iboot') ; max95 = bootaerm.quantile(0.975, dim='iboot')
min95_3 = bootaer3m.quantile(0.025, dim='iboot') ; max95_3 = bootaer3m.quantile(0.975, dim='iboot')
 
min95 = min95.rename('TOA_min95') ; max95 = max95.rename('TOA_max95')
min95_3 = min95_3.rename('TOA_min95_3') ; max95_3 = max95_3.rename('TOA_max95_3')

dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')

dat.to_netcdf(pathout+'CESM2_bootstrap_TOA.nc')
 
