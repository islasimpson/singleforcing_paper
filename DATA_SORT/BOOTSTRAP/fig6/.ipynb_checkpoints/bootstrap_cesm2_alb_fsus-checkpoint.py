import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial


from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig7/"

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds['time'] = ds.time.dt.year
    ds = ds.rename({'time':'year'})
    ds = ds.sel(year=slice(1920,2050))

    return ds

# 21 year running means
def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

aer2_fsns = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
           "AAER_FSNS_jja.nc", preprocess=partial(preprocessor)).FSNS
aer2_fsds = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
           "AAER_FSDS_jja.nc", preprocess=partial(preprocessor)).FSDS

aer_fsus = aer2_fsds - aer2_fsns
aer_fsus_gm = avg.cosweightlonlat(aer_fsus, 0, 360, 50, 90)
aer_fsds_gm = avg.cosweightlonlat(aer2_fsds, 0, 360, 50, 90)

aer_alb = aer_fsus_gm / aer_fsds_gm

aer_alb_base = aer_alb.sel(year=slice(1920,1940)).mean('year')
aer_fsus_base = aer_fsus_gm.sel(year=slice(1920,1940)).mean('year')

aer_alb = aer_alb - aer_alb_base
aer_fsus_gm = aer_fsus_gm - aer_fsus_base

aer_alb_21y = calc21ymean(aer_alb)
aer_fsus_21y = calc21ymean(aer_fsus_gm)

bootaeralb = boot.bootgen(aer_alb_21y, nsamples = aer_alb_21y.M.size, seed=3, nboots=1000)
bootaeralb3 = boot.bootgen(aer_alb_21y, nsamples=3, seed=4, nboots=1000)
bootaerfsus = boot.bootgen(aer_fsus_21y, nsamples = aer_fsus_21y.M.size, seed=3, nboots=1000)
bootaerfsus3 = boot.bootgen(aer_fsus_21y, nsamples=3, seed=4, nboots=1000)

bootaeralbm = bootaeralb.mean('isample')
bootaeralb3m = bootaeralb3.mean('isample')
bootaerfsusm = bootaerfsus.mean('isample')
bootaerfsus3m = bootaerfsus3.mean('isample')


min95alb = bootaeralbm.quantile(0.025, dim='iboot') ; max95alb = bootaeralbm.quantile(0.975, dim='iboot')
min95alb3 = bootaeralb3m.quantile(0.025, dim='iboot') ; max95alb3 = bootaeralb3m.quantile(0.975, dim='iboot')
min95fsus = bootaerfsusm.quantile(0.025, dim='iboot') ; max95fsus = bootaerfsusm.quantile(0.975, dim='iboot')
min95fsus3 = bootaerfsus3m.quantile(0.025, dim='iboot') ; max95fsus3 = bootaerfsus3m.quantile(0.975, dim='iboot')

min95alb = min95alb.rename('Albedo_min95') ; max95alb = max95alb.rename('Albedo_max95')
min95alb3 = min95alb3.rename('Albedo_min95_3') ; max95alb3 = max95alb3.rename('Albedo_max95_3')
min95fsus = min95fsus.rename('FSUS_min95') ; max95fsus = max95fsus.rename('FSUS_max95')
min95fsus3 = min95fsus3.rename('FSUS_min95_3') ; max95fsus3 = max95fsus3.rename('FSUS_max95_3')

dat = xr.merge([min95alb, max95alb, min95alb3, max95alb3,
               min95fsus, max95fsus, min95fsus3, max95fsus3], compat='override')

dat.to_netcdf(pathout+'CESM2_bootstrap_Albedo_FSUS.nc')

