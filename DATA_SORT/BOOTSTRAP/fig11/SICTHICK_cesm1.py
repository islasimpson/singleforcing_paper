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

dat = xr.open_dataset("/project/mojave/cesm2/LENS/ice/month_1/hi/b.e21.BSSP370smbb.f09_g17.LE2-1301.018.cice.h.hi.205501-206412.nc")
wgts = dat.tarea
wgts = xr.where( (dat.TLON > 220) & (dat.TLON < 360), wgts, 0)
wgts = xr.where( (dat.TLAT > 70) & (dat.TLAT < 90), wgts, 0)

lens1t = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/LENS1_sicthick_am.nc")
lens1t = lens1t.hi
lens1t = lens1t.transpose("M","year","nj","ni")
lens1 = xr.DataArray(np.zeros([lens1t.M.size, lens1t.year.size, wgts.nj.size, lens1t.ni.size]),
                     coords=[lens1t.M, lens1t.year, wgts.nj, lens1t.ni],
                     dims=['M','year','nj','ni'], name='hi')
lens1[:,:,384-lens1t.nj.size:384,:] = np.array(lens1t)

xaaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/XAER_sicthick_am.nc")
xaaer1 = xaaer1.hi

#------bootstrap the picontrol for AAER
picontrol = xr.open_dataset(
      "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM1/"+
      "SICTHICK_annualmean_piControl_northofgreenland.nc")
bootpi_base = boot.bootgenchunk_multimem(picontrol.hi, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(picontrol.hi, nyears=21, nmems=3, nboots=10000, seed=5)
bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])
#----------------------------------------

aaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"+
       "AAER1_sicthick_am.nc")
aaer1 = aaer1.hi


lens1_base = lens1.sel(year=slice(1920,1940)).mean('year')
xaaer1_base = xaaer1.sel(year=slice(1920,1940)).mean('year')
aaer1_base = aaer1.sel(year=slice(1920,1940)).mean('year')

lens1_anoms = lens1 - lens1_base
xaaer1_anoms = xaaer1 - xaaer1_base
aaer1_anoms = aaer1 - aaer1_base

lens1_anoms_w = lens1_anoms.weighted(wgts.fillna(0))
lens1_anoms_m = lens1_anoms_w.mean(("ni","nj"))

xaaer1_anoms_w = xaaer1_anoms.weighted(wgts.fillna(0))
xaaer1_anoms_m = xaaer1_anoms_w.mean(("ni","nj"))

aaer1_anoms_w = aaer1_anoms.weighted(wgts.fillna(0))
aaer1_anoms_m = aaer1_anoms_w.mean(("ni","nj"))

lens1_hi = lens1_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
aaer1_hi = aaer1_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
xaaer1_hi = xaaer1_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
xaaer1_hi = xaaer1_hi.transpose("M","year")

bootlens1 = boot.bootgen(lens1_hi, nsamples = lens1_hi.M.size, seed=4, nboots=1000)
bootlens1_3 = boot.bootgen(lens1_hi, nsamples = 3, seed=5, nboots=1000)

bootxaaer1 = boot.bootgen(xaaer1_hi, nsamples = xaaer1_hi.M.size, seed=4, nboots=1000)
bootxaaer1_3 = boot.bootgen(xaaer1_hi, nsamples = 3, seed=5, nboots=1000)

bootlens1m = bootlens1.mean('isample')
bootlens1_3m = bootlens1_3.mean('isample')
bootxaaer1m = bootxaaer1.mean('isample')
bootxaaer1_3m = bootxaaer1_3.mean('isample')


bootaaer1m = xr.DataArray(np.zeros([1000, aaer1_hi.year.size]),
                 coords=[np.arange(0,1000,1), aaer1_hi.year],
                 dims=['iboot','year'])
for iyear in np.arange(0,aaer1_hi.year.size,1):
    bootaaer1m[:,iyear] = aaer1_hi.mean('M').isel(year=iyear)+ (bootpi_period[:] - bootpi_base[:])

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
bootaer1xway3 = bootlens1m -  bootxaaer1_3m

min95_aer1xway = bootaer1xway.quantile(0.025, dim='iboot')
max95_aer1xway = bootaer1xway.quantile(0.975, dim='iboot')

min95_aer1xway3 = bootaer1xway3.quantile(0.025, dim='iboot')
max95_aer1xway3 = bootaer1xway3.quantile(0.975, dim='iboot')

min95_aer1xway = min95_aer1xway.rename('min95_aer1xway')
max95_aer1xway = max95_aer1xway.rename('max95_aer1xway')
min95_aer1xway3 = min95_aer1xway3.rename('min95_aer1xway3')
max95_aer1xway3 = max95_aer1xway3.rename('max95_aer1xway3')




dat = xr.merge([min95_lens1, max95_lens1, min95_lens1_3, max95_lens1_3,
               min95_aaer1, max95_aaer1, min95_xaer1, max95_xaer1,
               min95_aer1xway, max95_aer1xway, min95_aer1xway3, max95_aer1xway3], compat='override')

dat.to_netcdf(pathout+'CESM1_SICTHICK.nc')





