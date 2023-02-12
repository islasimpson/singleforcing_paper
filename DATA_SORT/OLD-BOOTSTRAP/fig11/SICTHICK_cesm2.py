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

lens2 = xr.open_dataset(
  "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_sicthick_am.nc")
lens2 = lens2.hi
aaer2 = xr.open_dataset(
  "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/AAER_sicthick_am.nc")
aaer2 = aaer2.hi
xaer2 = xr.open_dataset(
  "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/xAER_sicthick_am.nc")
xaer2 = xaer2.hi


#-----bootstrap the piControl for XAAER
picontrol = xr.open_dataset(
   "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/"+
   "SICTHICK_annualmean_picontrol_northofgreenland.nc")
bootpi_base = boot.bootgenchunk_multimem(picontrol.hi, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(picontrol.hi, nyears=21, nmems=3, nboots=1000, seed=5)
bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])
#----------------------------------------


lens2_base = lens2.sel(year=slice(1920,1940)).mean('year')
aaer2_base = aaer2.sel(year=slice(1920,1940)).mean('year')
xaer2_base = xaer2.sel(year=slice(1920,1940)).mean('year')

lens2_anoms = lens2 - lens2_base
aaer2_anoms = aaer2 - aaer2_base
xaer2_anoms = xaer2 - xaer2_base

lens2_anoms_w = lens2_anoms.weighted(wgts.fillna(0))
lens2_anoms_m = lens2_anoms_w.mean(("ni","nj"))

aaer2_anoms_w = aaer2_anoms.weighted(wgts.fillna(0))
aaer2_anoms_m = aaer2_anoms_w.mean(("ni","nj"))

xaer2_anoms_w = xaer2_anoms.weighted(wgts.fillna(0))
xaer2_anoms_m = xaer2_anoms_w.mean(("ni","nj"))

lens2_hi = lens2_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
aaer2_hi = aaer2_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
xaer2_hi = xaer2_anoms_m.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')

lens2_hi = lens2_hi.transpose("M","year")

bootlens2 = boot.bootgen(lens2_hi, nsamples = lens2_hi.M.size, seed=4, nboots=1000)
bootlens2_3 = boot.bootgen(lens2_hi, nsamples = 3, seed=5, nboots=1000)

bootaaer2 = boot.bootgen(aaer2_hi, nsamples = aaer2_hi.M.size, seed=4, nboots=1000)
bootaaer2_3 = boot.bootgen(aaer2_hi, nsamples = 3, seed=5, nboots=1000)

bootlens2m = bootlens2.mean('isample')
bootlens2_3m = bootlens2_3.mean('isample')
bootaaer2m = bootaaer2.mean('isample')
bootaaer2_3m = bootaaer2_3.mean('isample')

bootxaer2m = xr.DataArray(np.zeros([1000, xaer2_hi.year.size]),
                 coords=[np.arange(0,1000,1), xaer2_hi.year],
                 dims=['iboot','year'])
for iyear in np.arange(0,xaer2_hi.year.size,1):
   bootxaer2m[:,iyear] = \
     xaer2_hi.mean('M').isel(year=iyear) + \
     (bootpi_period[:] - bootpi_base[:])
       

min95_lens2 = bootlens2m.quantile(0.025, dim='iboot')
max95_lens2 = bootlens2m.quantile(0.975, dim='iboot')
min95_lens2_3 = bootlens2_3m.quantile(0.025, dim='iboot')
max95_lens2_3 = bootlens2_3m.quantile(0.975, dim='iboot')

min95_lens2 = min95_lens2.rename('min95_lens2')
max95_lens2 = max95_lens2.rename('max95_lens2')
min95_lens2_3 = min95_lens2_3.rename('min95_lens2_3')
max95_lens2_3 = max95_lens2_3.rename('max95_lens2_3')

min95_aaer2 = bootaaer2m.quantile(0.025, dim='iboot')
max95_aaer2 = bootaaer2m.quantile(0.975, dim='iboot')
min95_aaer2_3 = bootaaer2_3m.quantile(0.025, dim='iboot')
max95_aaer2_3 = bootaaer2_3m.quantile(0.975, dim='iboot')

min95_aaer2 = min95_aaer2.rename('min95_aaer2')
max95_aaer2 = max95_aaer2.rename('max95_aaer2')
min95_aaer2_3 = min95_aaer2_3.rename('min95_aaer2_3')
max95_aaer2_3 = max95_aaer2_3.rename('max95_aaer2_3')

min95_xaer2 = bootxaer2m.quantile(0.025, dim='iboot')
max95_xaer2 = bootxaer2m.quantile(0.975, dim='iboot')
min95_xaer2 = min95_xaer2.rename('min95_xaer2')
max95_xaer2 = max95_xaer2.rename('max95_xaer2')




dat = xr.merge([min95_lens2, max95_lens2, min95_lens2_3, max95_lens2_3,
               min95_aaer2, max95_aaer2, min95_aaer2_3, max95_aaer2_3,
               min95_xaer2, max95_xaer2], compat='override')

dat.to_netcdf(pathout+'CESM2_SICTHICK.nc')





