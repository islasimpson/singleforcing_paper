import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from CASutils import regrid_utils as regrid
from CASutils import mapplot_utils as mymaps
from CASutils import colorbar_utils as cbars
from CASutils import averaging_utils as avg
from CASutils import linfit_utils as linfit
from CASutils import bootstrap_utils as boot

import warnings
import sys
warnings.filterwarnings("ignore")

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig11/"


### Define atmosphere grid for regridding
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
atmgrid=xr.Dataset({'lat':landfrac.lat, 'lon':landfrac.lon})

# Read in the weights and set the values for the region outside the lab sea to zero
wgts=xr.open_dataset("/project/mojave/cesm2/LENS/ocn/month_1/HMXL/b.e21.BSSP370smbb.f09_g17.LE2-1301.018.pop.h.HMXL.208501-209412.nc").TAREA
wgts = xr.where( (wgts.TLONG > 300) & (wgts.TLONG < 315), wgts, 0)
wgts = xr.where( (wgts.TLAT > 53) & (wgts.TLAT < 65), wgts, 0)

lens1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/LENS1_March_HMXL.nc")
lens1['time'] = lens1.time.dt.year

xaaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/XAER_March_HMXL.nc")
xaaer1['time'] = xaaer1.time.dt.year

aaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"+
       "AAER_March_HMXL.nc")
aaer1 = aaer1.HMXL
aaer1['time'] = aaer1.time.dt.year

#------bootstrap the picontrol for AAER
picontrol = xr.open_dataset(
      "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM1/"+
      "HMXL_March_piControl_labsea.nc")
bootpi_base = boot.bootgenchunk_multimem(picontrol.HMXL, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(picontrol.HMXL, nyears=21, nmems=3, nboots=1000, seed=5)
bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])
#----------------------------------------


lens1base = lens1.sel(time=slice(1920,1940)).mean('time')
xaaer1base = xaaer1.sel(time=slice(1920,1940)).mean('time')
aaer1base = aaer1.sel(time=slice(1920,1940)).mean('time')

lens1_anoms = lens1 - lens1base
xaaer1_anoms = xaaer1 - xaaer1base
aaer1_anoms = aaer1 - aaer1base

lens1_anoms_w = lens1_anoms.weighted(wgts.fillna(0))
lens1_anoms_m = lens1_anoms_w.mean(("nlon","nlat"))

xaaer1_anoms_w = xaaer1_anoms.weighted(wgts.fillna(0))
xaaer1_anoms_m = xaaer1_anoms_w.mean(("nlon","nlat"))

aaer1_anoms_w = aaer1_anoms.weighted(wgts.fillna(0))
aaer1_anoms_m = aaer1_anoms_w.mean(("nlon","nlat"))

lens1_hmxl = lens1_anoms_m.rolling(time=21, min_periods=21, center='True').mean('time').dropna('time')
xaaer1_hmxl = xaaer1_anoms_m.rolling(time=21, min_periods=21, center='True').mean('time').dropna('time')
aaer1_hmxl = aaer1_anoms_m.rolling(time=21, min_periods=21, center='Trye').mean('time').dropna('time')

bootlens1 = boot.bootgen(lens1_hmxl.HMXL, nsamples = lens1_hmxl.M.size, seed=4, nboots=1000)
bootlens1_3 = boot.bootgen(lens1_hmxl.HMXL, nsamples = 3, seed=5, nboots=1000)

bootxaaer1 = boot.bootgen(xaaer1_hmxl.HMXL, nsamples = xaaer1_hmxl.M.size, seed=4, nboots=1000)
bootxaaer1_3 = boot.bootgen(xaaer1_hmxl.HMXL, nsamples = 3, seed=5, nboots=1000)

bootlens1m = bootlens1.mean('isample')
bootlens1_3m = bootlens1_3.mean('isample')
bootxaaer1m = bootxaaer1.mean('isample')
bootxaaer1_3m = bootxaaer1_3.mean('isample')

bootaaer1m = xr.DataArray(np.zeros([1000, aaer1_hmxl.time.size]),
                 coords=[np.arange(0,1000,1), aaer1_hmxl.time],
                 dims=['iboot','time'])
for iyear in np.arange(0,aaer1_hmxl.time.size,1):
    bootaaer1m[:,iyear] = aaer1_hmxl.mean('M').isel(time=iyear)+ (bootpi_period[:] - bootpi_base[:])




min95_lens1 = bootlens1m.quantile(0.025, dim='iboot')
max95_lens1 = bootlens1m.quantile(0.975, dim='iboot')
min95_lens1_3 = bootlens1_3m.quantile(0.025, dim='iboot')
max95_lens1_3 = bootlens1_3m.quantile(0.975, dim='iboot')

min95_lens1 = min95_lens1.rename('min95_lens1')
max95_lens1 = max95_lens1.rename('max95_lens1')
min95_lens1_3 = min95_lens1_3.rename('min95_lens1_3')
max95_lens1_3 = max95_lens1_3.rename('max95_lens1_3')

min95_aaer1 = bootaaer1m.quantile(0.025, dim='iboot')
max95_aaer1 = bootaaer1m.quantile(0.975, dim='iboot')

min95_aaer1 = min95_aaer1.rename('min95_aaer1')
max95_aaer1 = max95_aaer1.rename('max95_aaer1')

min95_xaer1 = bootxaaer1m.quantile(0.025, dim='iboot')
max95_xaer1 = bootxaaer1m.quantile(0.975, dim='iboot')

min95_xaer1 = min95_xaer1.rename('min95_xaer1')
max95_xaer1 = max95_xaer1.rename('max95_xaer1')

bootaer1xway = bootlens1m - bootxaaer1m
bootaer1xway3 = bootlens1m - bootxaaer1_3m

min95_aer1xway = bootaer1xway.quantile(0.025, dim='iboot')
max95_aer1xway = bootaer1xway.quantile(0.975, dim='iboot')

min95_aer1xway3 = bootaer1xway3.quantile(0.025, dim='iboot')
max95_aer1xway3 = bootaer1xway3.quantile(0.975, dim='iboot')

min95_aer1xway = min95_aer1xway.rename('min95_aer1xway')
max95_aer1xway = max95_aer1xway.rename('max95_aer1xway')
min95_aer1xway3 = min95_aer1xway3.rename('min95_aer1xway3')
max95_aer1xway3 = max95_aer1xway3.rename('max95_aer1xway3')



dat = xr.merge([min95_lens1, max95_lens1, min95_lens1_3, max95_lens1_3,
               min95_aaer1, max95_aaer1,min95_xaer1, max95_xaer1,
               min95_aer1xway, max95_aer1xway, min95_aer1xway3, max95_aer1xway3], compat='override')

dat.to_netcdf(pathout+'CESM1_HMXL.nc')





