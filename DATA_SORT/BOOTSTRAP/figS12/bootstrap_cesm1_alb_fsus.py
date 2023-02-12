import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial


from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/figS12/"

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

xaer1_fsns = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
       "LENS1-SF/XAER_FSNS_jja.nc", preprocess=partial(preprocessor)).FSNS
xaer1_fsds = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
       "LENS1-SF/XAER_FSDS_jja.nc", preprocess=partial(preprocessor)).FSDS
xaer1_fsus = xaer1_fsds - xaer1_fsns

lens1_fsns = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
       "LENS1/LENS1_FSNS_jja.nc", preprocess=partial(preprocessor)).FSNS
lens1_fsds = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
       "LENS1/LENS1_FSDS_jja.nc", preprocess=partial(preprocessor)).FSDS
lens1_fsus = lens1_fsds - lens1_fsns

xaer1_fsus_gm = avg.cosweightlonlat(xaer1_fsus, 0, 360, 50, 90)
xaer1_fsds_gm = avg.cosweightlonlat(xaer1_fsds, 0, 360, 50, 90)
xaer1_alb = xaer1_fsus_gm / xaer1_fsds_gm

lens1_fsus_gm = avg.cosweightlonlat(lens1_fsus, 0, 360, 50, 90)
lens1_fsds_gm = avg.cosweightlonlat(lens1_fsds, 0, 360, 50, 90)
lens1_alb = lens1_fsus_gm / lens1_fsds_gm

xaer1_alb_base = xaer1_alb.sel(year=slice(1920,1940)).mean('year')
xaer1_fsus_base = xaer1_fsus_gm.sel(year=slice(1920,1940)).mean('year')

lens1_alb_base = lens1_alb.sel(year=slice(1920,1940)).mean('year')
lens1_fsus_base = lens1_fsus_gm.sel(year=slice(1920,1940)).mean('year')

xaer1_alb = xaer1_alb - xaer1_alb_base
xaer1_fsus_gm = xaer1_fsus_gm - xaer1_fsus_base

lens1_alb = lens1_alb - lens1_alb_base
lens1_fsus_gm = lens1_fsus_gm - lens1_fsus_base

xaer1_alb_21y = calc21ymean(xaer1_alb)
xaer1_fsus_21y = calc21ymean(xaer1_fsus_gm)

lens1_alb_21y = calc21ymean(lens1_alb)
lens1_fsus_21y = calc21ymean(lens1_fsus_gm)

bootxaer1alb = boot.bootgen(xaer1_alb_21y, nsamples = xaer1_alb_21y.M.size, seed=3, nboots=1000)
bootxaer1alb3 = boot.bootgen(xaer1_alb_21y, nsamples = 3, seed=4, nboots=1000)
bootxaer1fsus = boot.bootgen(xaer1_fsus_21y, nsamples = xaer1_fsus_21y.M.size, seed=3, nboots=1000)
bootxaer1fsus3 = boot.bootgen(xaer1_fsus_21y, nsamples = 3, seed=4, nboots=1000)

bootlens1alb = boot.bootgen(lens1_alb_21y, nsamples = lens1_alb_21y.M.size, seed=3, nboots=1000)
bootlens1fsus = boot.bootgen(lens1_fsus_21y, nsamples = lens1_fsus_21y.M.size, seed=3, nboots=1000)


bootxaer1albm = bootxaer1alb.mean('isample')
bootxaer1alb3m = bootxaer1alb3.mean('isample')
bootxaer1fsusm = bootxaer1fsus.mean('isample')
bootxaer1fsus3m = bootxaer1fsus3.mean('isample')
bootlens1albm = bootlens1alb.mean('isample')
bootlens1fsusm = bootlens1fsus.mean('isample')

bootaer1xwayalb = bootlens1albm - bootxaer1albm
bootaer1xwayalb3 = bootlens1albm - bootxaer1alb3m
bootaer1xwayfsus = bootlens1fsusm - bootxaer1fsusm
bootaer1xwayfsus3 = bootlens1fsusm - bootxaer1fsus3m

min95alb = bootaer1xwayalb.quantile(0.025, dim='iboot')
max95alb = bootaer1xwayalb.quantile(0.975, dim='iboot')
min95alb3 = bootaer1xwayalb3.quantile(0.025, dim='iboot')
max95alb3 = bootaer1xwayalb3.quantile(0.975, dim='iboot')

min95fsus = bootaer1xwayfsus.quantile(0.025, dim='iboot')
max95fsus = bootaer1xwayfsus.quantile(0.975, dim='iboot')
min95fsus3 = bootaer1xwayfsus3.quantile(0.025, dim='iboot')
max95fsus3 = bootaer1xwayfsus3.quantile(0.975, dim='iboot')

min95alb = min95alb.rename('Albedo_min95') ; max95alb = max95alb.rename('Albedo_max95')
min95alb3 = min95alb3.rename('Albedo_min95_3') ; max95alb3 = max95alb3.rename('Albedo_max95_3')
min95fsus = min95fsus.rename('FSUS_min95') ; max95fsus = max95fsus.rename('FSUS_max95')
min95fsus3 = min95fsus3.rename('FSUS_min95_3') ; max95fsus3 = max95fsus3.rename('FSUS_max95_3')

dat = xr.merge([min95alb, max95alb, min95alb3, max95alb3,
               min95fsus, max95fsus, min95fsus3, max95fsus3], compat='override')

dat.to_netcdf(pathout+'CESM1_bootstrap_Albedo_FSUS.nc')











