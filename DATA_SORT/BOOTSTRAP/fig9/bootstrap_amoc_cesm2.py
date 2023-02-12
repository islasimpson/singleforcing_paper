import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

from math import nan

from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit


pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig9/"

def calc21ymean(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

#-------------LENS2
all2_amoc = xr.open_dataset(
     "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_AMOC45_am.nc")
all2_amoc = all2_amoc.MOC
all2_amoc = all2_amoc.transpose("M","year")
all2_base = all2_amoc.sel(year=slice(1920,1940)).mean('year')
all2anoms = all2_amoc - all2_base
all2anoms_sm = calc21ymean(all2anoms)
#------------------

#-------------AAER
aaer2_amoc = xr.open_dataset(
     "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/AAER_AMOC45_am.nc")
aaer2_amoc = aaer2_amoc.MOC
aaer2_amoc = aaer2_amoc.transpose("M","year")

aaer2_base = aaer2_amoc.sel(year=slice(1920,1940)).mean('year')
aaer2anoms = aaer2_amoc - aaer2_base
aaer2anoms_sm = calc21ymean(aaer2anoms)
#-----------------

#-------------GHG
ghg2_amoc = xr.open_dataset(
     "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/GHG_AMOC45_am.nc")
ghg2_amoc = ghg2_amoc.MOC
ghg2_amoc = ghg2_amoc.transpose("M","year")

ghg2_base = ghg2_amoc.sel(year=slice(1920,1940)).mean('year')
ghg2anoms = ghg2_amoc - ghg2_base
ghg2anoms_sm = calc21ymean(ghg2anoms)
#-----------------

#-------------BMB
bmb2_amoc = xr.open_dataset(
     "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/BMB_AMOC45_am.nc")
bmb2_amoc = bmb2_amoc.MOC
bmb2_amoc = bmb2_amoc.transpose("M","year")

bmb2_base = bmb2_amoc.sel(year=slice(1920,1940)).mean('year')
bmb2anoms = bmb2_amoc - bmb2_base
bmb2anoms_sm = calc21ymean(bmb2anoms)
#-----------------

#-------------EE
ee2_amoc = xr.open_dataset(
     "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/EE_AMOC45_am.nc")
ee2_amoc = ee2_amoc.MOC
ee2_amoc = ee2_amoc.transpose("M","year")

ee2_base = ee2_amoc.sel(year=slice(1920,1940)).mean('year')
ee2anoms = ee2_amoc - ee2_base
ee2anoms_sm = calc21ymean(ee2anoms)
#-----------------

#-------------XAER
xaer2_amoc = xr.open_dataset(
    "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/xAER_AMOC45_am.nc")
xaer2_amoc = xaer2_amoc.MOC
xaer2_amoc = xaer2_amoc.transpose("M","year")
xaer2_base = xaer2_amoc.sel(year=slice(1920,1940)).mean('year')
xaer2anoms = xaer2_amoc - xaer2_base
xaer2anoms_sm = calc21ymean(xaer2anoms)
#-----------------


#-------piControl
picontrol = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/"+
   "piControl/CESM2/piControl_AMOC_am.nc")
picontrol = picontrol.MOC
a, b = linfit.linfit_xy(picontrol.year, picontrol)
picontrol = picontrol - (a + b*picontrol.year)
picontrol = calc21ymean(picontrol)


bootall = boot.bootgen(all2anoms_sm, nsamples = all2_amoc.M.size, seed=4, nboots=1000)
bootghg = boot.bootgen(ghg2anoms_sm, nsamples = ghg2_amoc.M.size, seed=4, nboots=1000)
bootaer = boot.bootgen(aaer2anoms_sm, nsamples = aaer2_amoc.M.size, seed=4, nboots=1000)
bootbmb = boot.bootgen(bmb2anoms_sm, nsamples = bmb2_amoc.M.size, seed=4, nboots=1000)
bootee = boot.bootgen(ee2anoms_sm, nsamples = ee2_amoc.M.size, seed=4, nboots=1000)
bootaer3 = boot.bootgen(aaer2anoms_sm, nsamples = 3, seed=5, nboots=1000)

#-----bootstrap the piControl to get a uncertainty on the XAAER simulation
bootxaaerbase = boot.bootgenchunk_multimem(picontrol, 1, 3, nboots=1000, seed=4)
bootxaaeranoms = boot.bootgenchunk_multimem(picontrol, 1, 3, nboots=1000, seed=5)

bootxaaerbasem = bootxaaerbase.mean('imem')
bootxaaeranomsm = bootxaaeranoms.mean('imem')
bootxaer = bootxaaeranomsm.isel(isample=0) - bootxaaerbasem.isel(isample=0)

bootallm = bootall.mean('isample')
bootghgm = bootghg.mean('isample')
bootaerm = bootaer.mean('isample')
bootbmbm = bootbmb.mean('isample')
booteem = bootee.mean('isample')
bootaer3m = bootaer3.mean('isample')

bootxaerm = xr.DataArray(np.zeros([1000, xaer2anoms_sm.year.size]),
               coords=[np.arange(0,1000,1), xaer2anoms_sm.year],
               dims=['iboot','year'])
for iyear in np.arange(0,xaer2anoms_sm.year.size,1):
    bootxaerm[:,iyear] = xaer2anoms_sm.isel(year=iyear).mean('M').values + bootxaer


bootghgplusaerm = bootghgm + bootaerm

min95_all = bootallm.quantile(0.025, dim='iboot')
max95_all = bootallm.quantile(0.975, dim='iboot')

min95_ghg = bootghgm.quantile(0.025, dim='iboot')
max95_ghg = bootghgm.quantile(0.975, dim='iboot')

min95_aaer = bootaerm.quantile(0.025, dim='iboot') 
max95_aaer = bootaerm.quantile(0.975, dim='iboot')

min95_bmb = bootbmbm.quantile(0.025, dim='iboot') 
max95_bmb = bootbmbm.quantile(0.975, dim='iboot')

min95_ee = booteem.quantile(0.025, dim='iboot') 
max95_ee = booteem.quantile(0.975, dim='iboot')

min95_aaer3 = bootaer3m.quantile(0.025, dim='iboot')
max95_aaer3 = bootaer3m.quantile(0.975, dim='iboot')

min95_xaer = bootxaerm.quantile(0.025, dim='iboot')
max95_xaer = bootxaerm.quantile(0.975, dim='iboot')

min95_ghgplusaer = bootghgplusaerm.quantile(0.025, dim='iboot')
max95_ghgplusaer = bootghgplusaerm.quantile(0.975, dim='iboot')

min95_all = min95_all.rename('min95_all')
max95_all = max95_all.rename('max95_all')
min95_ghg = min95_ghg.rename('min95_ghg')
max95_ghg = max95_ghg.rename('max95_ghg')
min95_aaer = min95_aaer.rename('min95_aaer')
max95_aaer = max95_aaer.rename('max95_aaer')
min95_bmb = min95_bmb.rename('min95_bmb')
max95_bmb = max95_bmb.rename('max95_bmb')
min95_ee = min95_ee.rename('min95_ee')
max95_ee = max95_ee.rename('max95_ee')
min95_aaer3 = min95_aaer3.rename('min95_aaer3')
max95_aaer3 = max95_aaer3.rename('max95_aaer3')
min95_xaer = min95_xaer.rename('min95_xaer')
max95_xaer = max95_xaer.rename('max95_xaer')
min95_ghgplusaer = min95_ghgplusaer.rename('min95_ghgplusaer')
max95_ghgplusaer = max95_ghgplusaer.rename('max95_ghgplusaer')


dat = xr.merge([min95_all, max95_all, min95_ghg, max95_ghg, min95_aaer, max95_aaer, 
       min95_bmb, max95_bmb, min95_ee, max95_ee,min95_ghgplusaer, max95_ghgplusaer,
       min95_aaer3, max95_aaer3, min95_xaer, max95_xaer], compat='override')

dat.to_netcdf(pathout+'AMOC_bootstrap_cesm2.nc')

