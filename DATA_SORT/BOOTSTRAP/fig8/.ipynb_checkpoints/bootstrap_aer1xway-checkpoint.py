import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from functools import partial

import dask
from CASutils import bootstrap_utils as boot

dask.config.set(**{'array.slicing.split_large_chunks': False})

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig9/"

landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
# pre-processor to ensure all lon and lat coordinates are the same.  Also adding flexibility for reading in seasonal data
# Taking the ensemble mean
def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds = ds.sel(year=slice(1920,2050))
    return ds

seas = 'am' ; var=['TREFHT']

baselens1 = '/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/'
lens1 = []
for ivar in var:
    dat = xr.open_mfdataset(baselens1+'LENS1_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    lens1.append(dat)
lens1 = xr.merge(lens1)

basexaer1 = '/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/'
xaer1 = []
for ivar in var:
    dat = xr.open_mfdataset(basexaer1+'XAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    xaer1.append(dat)    
xaer1 = xr.merge(xaer1)

lens1base = lens1.sel(year=slice(1920,1940)).mean('year')
xaer1base = xaer1.sel(year=slice(1920,1940)).mean('year')

lens1period = lens1.sel(year=slice(2030,2050)).mean('year')
xaer1period = xaer1.sel(year=slice(2030,2050)).mean('year')


lens1anoms = lens1period - lens1base
xaer1anoms = xaer1period - xaer1base

bootslens1 = boot.bootgen(lens1anoms.TREFHT, nsamples = lens1anoms.M.size, seed=3, nboots=1000)
bootsxaer1 = boot.bootgen(xaer1anoms.TREFHT, nsamples = xaer1anoms.M.size, seed=4, nboots=1000)

bootslens1m = bootslens1.mean('isample')
bootsxaer1m = bootsxaer1.mean('isample')

bootsaer1xway = bootslens1m - bootsxaer1m

min95 = bootsaer1xway.quantile(0.025, dim='iboot')
max95 = bootsaer1xway.quantile(0.975, dim='iboot')

min95 = min95.rename('min95')
max95 = max95.rename('max95')

min95.to_netcdf(pathout+'aaer1xway_ci.nc')
max95.to_netcdf(pathout+'aaer1xway_ci.nc', mode='a')
