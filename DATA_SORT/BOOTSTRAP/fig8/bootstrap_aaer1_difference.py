import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from functools import partial

import dask
from CASutils import bootstrap_utils as boot

dask.config.set(**{'array.slicing.split_large_chunks': False})

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig8/"

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

baseaaer1 = '/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/'
aaer1 = []
for ivar in var:
    dat = xr.open_mfdataset(baseaaer1+'AAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    aaer1.append(dat)
aaer1 = xr.merge(aaer1)

aaer1base = aaer1.sel(year=slice(1920,1940)).mean('year')
lens1base = lens1.sel(year=slice(1920,1940)).mean('year')
xaer1base = xaer1.sel(year=slice(1920,1940)).mean('year')

aaer1period = aaer1.sel(year=slice(2030,2050)).mean('year')
lens1period = lens1.sel(year=slice(2030,2050)).mean('year')
xaer1period = xaer1.sel(year=slice(2030,2050)).mean('year')

lens1anoms = lens1period - lens1base
xaer1anoms = xaer1period - xaer1base
aaer1anoms = aaer1period - aaer1base

bootlens1 = boot.bootgen(lens1anoms.TREFHT, nsamples = lens1anoms.M.size, seed=3, nboots=1000)
bootxaer1 = boot.bootgen(xaer1anoms.TREFHT, nsamples = xaer1anoms.M.size, seed=4, nboots=1000)

bootlens1m = bootlens1.mean('isample')
bootxaer1m = bootxaer1.mean('isample')

baseaaer1 = '/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/'
aaer1 = []
for ivar in var:
    dat = xr.open_mfdataset(baseaaer1+'AAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    aaer1.append(dat)
aaer1 = xr.merge(aaer1)

basepi = '/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM1/'
pi = xr.open_dataset(basepi+'picontrol_TREFHT_am.nc')
pi = pi.TREFHT
pi['lon'] = landfrac.lon ; pi['lat'] = landfrac.lat

bootpi_base = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=5)

bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])

bootaaer1 = aaer1anoms.TREFHT.mean('M') + (bootpi_period - bootpi_base)
bootaaer1xway = bootlens1m - bootxaer1m

bootdiff = bootaaer1 - bootaaer1xway

min95 = bootdiff.quantile(0.025, dim='iboot')
max95 = bootdiff.quantile(0.975, dim='iboot')

min95 = min95.rename('min95')
max95 = max95.rename('max95')

min95.to_netcdf(pathout+'aer1diff_ci.nc')
max95.to_netcdf(pathout+'aer1diff_ci.nc', mode='a')


