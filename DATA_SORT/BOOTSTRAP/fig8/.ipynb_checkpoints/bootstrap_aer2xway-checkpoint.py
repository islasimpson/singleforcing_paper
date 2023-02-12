import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from functools import partial

import sys
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

#------------------CESM2
baselens2='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/'
lens2 = []
for ivar in var:
    dat = xr.open_mfdataset(baselens2+'LENS2_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    lens2.append(dat)
lens2 = xr.merge(lens2)

#------------------XAER2
basexaer2="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"
xaer2 = []
for ivar in var:
    dat = xr.open_mfdataset(basexaer2+'xAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    xaer2.append(dat)
xaer2 = xr.merge(xaer2)

#------------------piControl
basepi='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/'
pi = xr.open_dataset(basepi+'picontrol_TREFHT_am.nc')
pi = pi.TREFHT


lens2base = lens2.sel(year=slice(1920,1940)).mean('year')
lens2period = lens2.sel(year=slice(2030,2050)).mean('year')
xaer2base = xaer2.sel(year=slice(1920,1940)).mean('year')
xaer2period = xaer2.sel(year=slice(2030,2050)).mean('year')

lens2anoms = lens2period - lens2base
xaer2anoms = xaer2period - xaer2base

bootlens2 = boot.bootgen(lens2anoms.TREFHT, nsamples=lens2anoms.M.size, seed=3, nboots=1000)

#---bootstrapping the picontrol to get bootstrapped anomalies for XAER2
bootpi_base = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=5)

bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])

bootlens2 = bootlens2.mean('isample')

xaer2anoms['lon'] = bootpi_base.lon
xaer2anoms['lat'] = bootpi_base.lat
bootxaer2 = xaer2anoms.TREFHT.mean('M') + (bootpi_period - bootpi_base)

bootlens2['lon'] = bootpi_base.lon
bootlens2['lat'] = bootpi_base.lat

bootaer2xway = bootlens2 - bootxaer2

min95 = bootaer2xway.quantile(0.025, dim='iboot')
max95 = bootaer2xway.quantile(0.975, dim='iboot')

min95 = min95.rename('min95')
max95 = max95.rename('max95')

min95.to_netcdf(pathout+'aer2xway_ci.nc')
max95.to_netcdf(pathout+'aer2xway_ci.nc', mode='a')
