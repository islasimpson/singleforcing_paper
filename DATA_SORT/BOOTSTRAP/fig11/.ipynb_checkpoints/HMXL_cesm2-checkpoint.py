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

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig12/"


### Define atmosphere grid for regridding
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
atmgrid=xr.Dataset({'lat':landfrac.lat, 'lon':landfrac.lon})

# Read in the weights and set the values for the region outside the lab sea to zero
wgts=xr.open_dataset("/project/mojave/cesm2/LENS/ocn/month_1/HMXL/b.e21.BSSP370smbb.f09_g17.LE2-1301.018.pop.h.HMXL.208501-209412.nc").TAREA
wgts = xr.where( (wgts.TLONG > 300) & (wgts.TLONG < 315), wgts, 0)
wgts = xr.where( (wgts.TLAT > 53) & (wgts.TLAT < 65), wgts, 0)

lens2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_March_HMXL.nc")
lens2['time'] = lens2.time.dt.year

aaer2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/AAER_March_HMXL.nc")
aaer2['time'] = aaer2.time.dt.year

xaer2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"+
      "xAER_March_HMXL.nc")
xaer2['time'] = xaer2.time.dt.year

#-----bootstrap the piControl for XAAER
picontrol = xr.open_dataset(
   "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/"+
   "HMXL_March_piControl_labsea.nc")
bootpi_base = boot.bootgenchunk_multimem(picontrol.HMXL, nyears=21, nmems=3, nboots=1000, seed=4)
bootpi_period = boot.bootgenchunk_multimem(picontrol.HMXL, nyears=21, nmems=3, nboots=1000, seed=5)
bootpi_base = bootpi_base.mean(['imem','isample'])
bootpi_period = bootpi_period.mean(['imem','isample'])
#----------------------------------------




lens2base = lens2.sel(time=slice(1920,1940)).mean('time')
aaer2base = aaer2.sel(time=slice(1920,1940)).mean('time')
xaer2base = xaer2.sel(time=slice(1920,1940)).mean('time')

lens2_anoms = lens2 - lens2base
aaer2_anoms = aaer2 - aaer2base
xaer2_anoms = xaer2 - xaer2base

lens2_anoms_w = lens2_anoms.weighted(wgts.fillna(0))
lens2_anoms_m = lens2_anoms_w.mean(("nlon","nlat"))

aaer2_anoms_w = aaer2_anoms.weighted(wgts.fillna(0))
aaer2_anoms_m = aaer2_anoms_w.mean(("nlon","nlat"))

xaer2_anoms_w = xaer2_anoms.weighted(wgts.fillna(0))
xaer2_anoms_m = xaer2_anoms_w.mean(("nlon","nlat"))

lens2_hmxl = lens2_anoms_m.rolling(time=21, min_periods=21, center='True').mean('time').dropna('time')
aaer2_hmxl = aaer2_anoms_m.rolling(time=21, min_periods=21, center='True').mean('time').dropna('time')
xaer2_hmxl = xaer2_anoms_m.rolling(time=21, min_periods=21, center='True').mean('time').dropna('time')

bootlens2 = boot.bootgen(lens2_hmxl.HMXL, nsamples = lens2_hmxl.M.size, seed=4, nboots=1000)
bootlens2_3 = boot.bootgen(lens2_hmxl.HMXL, nsamples = 3, seed=5, nboots=1000)

bootaaer2 = boot.bootgen(aaer2_hmxl.HMXL, nsamples = aaer2_hmxl.M.size, seed=4, nboots=1000)
bootaaer2_3 = boot.bootgen(aaer2_hmxl.HMXL, nsamples = 3, seed=5, nboots=1000)

bootlens2m = bootlens2.mean('isample')
bootlens2_3m = bootlens2_3.mean('isample')
bootaaer2m = bootaaer2.mean('isample')
bootaaer2_3m = bootaaer2_3.mean('isample')

bootxaer2m = xr.DataArray(np.zeros([1000, xaer2_hmxl.time.size]),
                 coords=[np.arange(0,1000,1), xaer2_hmxl.time],
                 dims=['iboot','time'])
for iyear in np.arange(0,xaer2_hmxl.time.size,1):
   bootxaer2m[:,iyear] = \
     xaer2_hmxl.HMXL.mean('M').isel(time=iyear) + \
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

dat.to_netcdf(pathout+'CESM2_HMXL.nc')





