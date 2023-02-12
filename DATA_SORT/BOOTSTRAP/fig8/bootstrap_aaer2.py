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

baselens2='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/'
lens2 = []
for ivar in var:
    dat = xr.open_mfdataset(baselens2+'AAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    lens2.append(dat)
lens2 = xr.merge(lens2)

lens2base = lens2.sel(year=slice(1920,1940)).mean('year')
lens2period = lens2.sel(year=slice(2030,2050)).mean('year')

lens2anoms = lens2period - lens2base

boots = boot.bootgen(lens2anoms.TREFHT, nsamples=lens2anoms.M.size, seed=3, nboots=1000)

bootsm = boots.mean('isample')

min95 = bootsm.quantile(0.025, dim='iboot')
max95 = bootsm.quantile(0.975, dim='iboot')

min95 = min95.rename('min95')
max95 = max95.rename('max95')

min95.to_netcdf(pathout+'aaer2_ci.nc')
max95.to_netcdf(pathout+'aaer2_ci.nc', mode='a')
