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

aaer1base = aaer1.sel(year=slice(1920,1940)).mean('year')
aaer1period = aaer1.sel(year=slice(2030,2050)).mean('year')
aaer1anoms = aaer1period - aaer1base

bootpi_base = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(pi, nyears=21, nmems=3, nboots=1000, seed=5)

bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])

bootaaer1 = aaer1anoms.TREFHT.mean('M') + (bootpi_period - bootpi_base)

min95 = bootaaer1.quantile(0.025, dim='iboot')
max95 = bootaaer1.quantile(0.975, dim='iboot')

min95 = min95.rename('min95')
max95 = max95.rename('max95')

min95.to_netcdf(pathout+'aer1_ci.nc')
max95.to_netcdf(pathout+'aer1_ci.nc', mode='a')


