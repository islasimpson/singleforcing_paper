import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

from math import nan

from CASutils import linfit_utils as linfit
from CASutils import bootstrap_utils as boot

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig9/"

def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

all1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+
         "LENS1_AMOC45_am.nc")
all1 = all1.MOC
all1 = all1.transpose("M","year")

xaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+
    "XAER_AMOC45_am.nc")
xaer1 = xaer1.MOC
xaer1 = xaer1.transpose("M","year")

xghg1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+
    "XGHG_AMOC45_am.nc")
xghg1 = xghg1.MOC
xghg1 = xghg1.transpose("M","year")

xbmb1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+
    "XBMB_AMOC45_am.nc")
xbmb1 = xbmb1.MOC
xbmb1 = xbmb1.transpose("M","year")

aaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM1-AAER/"+
   "AAER_amoc_45_am.nc")
aaer1 = aaer1.MOC
aaer1 = aaer1.transpose("M","year")

picontrol = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
            "piControl/CESM1/piControl_AMOC_am.nc")
picontrol = picontrol.MOC
a, b = linfit.linfit_xy(picontrol.year, picontrol)
picontrol = picontrol - (a + b*picontrol.year)



all1base = all1.sel(year=slice(1920,1940)).mean('year')
xaer1base = xaer1.sel(year=slice(1920,1940)).mean('year')
xghg1base = xghg1.sel(year=slice(1920,1940)).mean('year')
xbmb1base = xbmb1.sel(year=slice(1920,1940)).mean('year')
aaer1base = aaer1.sel(year=slice(1920,1940)).mean('year')

all1 = all1 - all1base
xaer1 = xaer1 - xaer1base
xghg1 = xghg1 - xghg1base
xbmb1 = xbmb1 - xbmb1base
aaer1 = aaer1 - aaer1base

all1 = calc21ymean(all1)
xaer1 = calc21ymean(xaer1)
xghg1 = calc21ymean(xghg1)
xbmb1 = calc21ymean(xbmb1)
aaer1 = calc21ymean(aaer1)
picontrol = calc21ymean(picontrol)

bootall = boot.bootgen(all1, nsamples = all1.M.size, seed=3, nboots=1000)
bootxghg = boot.bootgen(xghg1, nsamples = xghg1.M.size, seed=4, nboots=1000)
bootxaer = boot.bootgen(xaer1, nsamples = xaer1.M.size, seed=4, nboots=1000)
bootxaer3 = boot.bootgen(xaer1, nsamples = 3, seed=4, nboots=1000)
bootxbmb = boot.bootgen(xbmb1, nsamples = xbmb1.M.size, seed=4, nboots=1000)

#-----bootstrap the piControl to get a ci on AAER
bootaaerbase = boot.bootgenchunk_multimem(picontrol, 1, 3, nboots=1000, seed=4)
bootaaeranoms = boot.bootgenchunk_multimem(picontrol, 1, 3, nboots=1000, seed=5)

bootaaerbasem = bootaaerbase.mean('imem')
bootaaeranomsm = bootaaeranoms.mean('imem')
bootaaer = bootaaeranomsm.isel(isample=0) - bootaaerbasem.isel(isample=0)

min95aaert = bootaaer.quantile(0.025, dim='iboot')
max95aaert = bootaaer.quantile(0.975, dim='iboot')

min95aaer = xr.DataArray(np.zeros([aaer1.year.size]), coords=[aaer1.year], dims=['year'], 
               name='min95_aaer')
max95aaer = xr.DataArray(np.zeros([aaer1.year.size]), coords=[aaer1.year], dims=['year'],
               name='max95_aaer')

min95aaer[:] = np.array(min95aaert)+aaer1.mean('M')
max95aaer[:] = np.array(max95aaert)+aaer1.mean('M')
#-------------------------------------

bootallm = bootall.mean('isample')
bootxghgm = bootxghg.mean('isample')
bootxaerm = bootxaer.mean('isample')
bootxaer3m = bootxaer3.mean('isample')
bootxbmbm = bootxbmb.mean('isample')

bootghgmxway = bootallm - bootxghgm
bootaermxway = bootallm - bootxaerm 
bootaermxway3 = bootallm - bootxaer3m
bootbmbmxway = bootallm - bootxbmbm

bootee = bootallm - (bootghgmxway + bootaermxway + bootbmbmxway)

bootghgplusaer = bootghgmxway + bootaermxway

min95_all = bootallm.quantile(0.025, dim='iboot') ; max95_all = bootallm.quantile(0.975, dim='iboot')
min95_xaer = bootxaerm.quantile(0.025, dim='iboot') ; max95_xaer = bootxaerm.quantile(0.975, dim='iboot')
min95_xghg = bootxghgm.quantile(0.025, dim='iboot') ; max95_xghg = bootxghgm.quantile(0.975, dim='iboot')
min95_xbmb = bootxbmbm.quantile(0.025, dim='iboot') ; max95_xbmb = bootxbmbm.quantile(0.975, dim='iboot')
min95_ee = bootee.quantile(0.025, dim='iboot') ; max95_ee = bootee.quantile(0.975, dim='iboot')

min95_aaerxway = bootaermxway.quantile(0.025, dim='iboot') 
max95_aaerxway = bootaermxway.quantile(0.975, dim='iboot')
min95_aaerxway3 = bootaermxway3.quantile(0.025, dim='iboot')
max95_aaerxway3 = bootaermxway3.quantile(0.975, dim='iboot')
min95_ghgxway = bootghgmxway.quantile(0.025, dim='iboot') 
max95_ghgxway = bootghgmxway.quantile(0.975, dim='iboot')
min95_bmbxway = bootbmbmxway.quantile(0.025, dim='iboot')
max95_bmbxway = bootbmbmxway.quantile(0.975, dim='iboot')

min95_ghgplusaer = bootghgplusaer.quantile(0.025, dim='iboot')
max95_ghgplusaer = bootghgplusaer.quantile(0.975, dim='iboot')


min95_all = min95_all.rename('min95_all')
max95_all = max95_all.rename('max95_all')
min95_xaer = min95_xaer.rename('min95_xaer')
max95_xaer = max95_xaer.rename('max95_xaer')
min95_xghg = min95_xghg.rename('min95_xghg')
max95_xghg = max95_xghg.rename('max95_xghg')
min95_ee = min95_ee.rename('min95_ee')
max95_ee = max95_ee.rename('max95_ee')


min95_aaerxway = min95_aaerxway.rename('min95_aaerxway')
max95_aaerxway = max95_aaerxway.rename('max95_aaerxway')
min95_aaerxway3 = min95_aaerxway3.rename('min95_aaerxway3')
max95_aaerxway3 = max95_aaerxway3.rename('max95_aaerxway3')
min95_ghgxway = min95_ghgxway.rename('min95_ghgxway')
max95_ghgxway = max95_ghgxway.rename('max95_ghgxway')
min95_bmbxway = min95_bmbxway.rename('min95_bmbxway')
max95_bmbxway = max95_bmbxway.rename('max95_bmbxway')


min95_ghgplusaer = min95_ghgplusaer.rename('min95_ghgplusaer')
max95_ghgplusaer = max95_ghgplusaer.rename('max95_ghgplusaer')

dat = xr.merge([min95_all, max95_all, min95_xaer, max95_xaer, min95_xghg, max95_xghg, 
         min95_ghgxway, max95_ghgxway, min95_aaerxway3, max95_aaerxway3, 
         min95aaer, max95aaer, min95_aaerxway, max95_aaerxway, 
         min95_bmbxway, max95_bmbxway, min95_ee, max95_ee,
         min95_ghgplusaer, max95_ghgplusaer], compat='override')

dat.to_netcdf(pathout+'AMOC_bootstrap_cesm1.nc')

