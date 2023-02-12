import xarray as xr
import numpy as np
import sys
from CASutils import linfit_utils as linfit
from CASutils import bootstrap_utils as boot

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig10/"

pi_amoc = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/piControl_AMOC_am.nc").MOC
pi_fsus = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/piControl_fsus_integral_am.nc").FSUS

#------detrend
a, b = linfit.linfit_xy(pi_fsus.year, pi_fsus)
pi_fsus = pi_fsus - (a + b*pi_fsus.year)

a, b = linfit.linfit_xy(pi_amoc.year, pi_amoc)
pi_amoc = pi_amoc - (a + b*pi_amoc.year)

pi_amoc_sm = pi_amoc.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
pi_fsus_sm = pi_fsus.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')


# Experiment AMOC's
aaer2_amoc = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/AAER_AMOC45_am.nc")
lens2_amoc = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_AMOC45_am.nc")
xaer2_amoc = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/xAER_AMOC45_am.nc")

aaer2_amoc_base = aaer2_amoc.sel(year=slice(1920,1940)).mean('year')
lens2_amoc_base = lens2_amoc.sel(year=slice(1920,1940)).mean('year')
xaer2_amoc_base = xaer2_amoc.sel(year=slice(1920,1940)).mean('year')

aaer2_amoc = aaer2_amoc - aaer2_amoc_base
lens2_amoc = lens2_amoc - lens2_amoc_base
xaer2_amoc = xaer2_amoc - xaer2_amoc_base

aaer2_amoc = aaer2_amoc.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
lens2_amoc = lens2_amoc.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
xaer2_amoc = xaer2_amoc.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')

aaer2boot = boot.bootgen(aaer2_amoc.MOC.transpose('M','year'), nsamples = aaer2_amoc.M.size, seed=4, nboots=1000)
lens2boot = boot.bootgen(lens2_amoc.MOC.transpose('M','year'), nsamples = lens2_amoc.M.size, seed=5, nboots=1000)

bootxaaerbase = boot.bootgenchunk_multimem(pi_amoc_sm, 1, 3, nboots=1000, seed=4)
bootxaaeranoms = boot.bootgenchunk_multimem(pi_amoc_sm, 1, 3, nboots=1000, seed=5)

bootxaaerbasem = bootxaaerbase.mean('imem')
bootxaaeranomsm = bootxaaeranoms.mean('imem')
bootxaer = bootxaaeranomsm.isel(isample=0) - bootxaaerbasem.isel(isample=0)


xaer2boot = xr.DataArray(np.zeros([1000, xaer2_amoc.year.size]),
                   coords=[np.arange(0,1000,1), xaer2_amoc.year],
                   dims=['iboot','year'])
for iyear in np.arange(0,xaer2_amoc.year.size, 1):
    xaer2boot[:,iyear] = xaer2_amoc.MOC.isel(year=iyear).mean('M').values + bootxaer 

aaer2bootm = aaer2boot.mean('isample')
lens2bootm = lens2boot.mean('isample')
#xaer2bootm = xaer2boot.mean('isample')

aaer2xway = lens2bootm - xaer2boot

amocdif = aaer2bootm - aaer2xway



amocboot = boot.bootgenchunk_multimem(pi_amoc_sm, 200, 10, 1000, seed=3)
amocboot = amocboot.stack(year=("imem","isample"))
amocboot['year'] = np.arange(0,2000,1)

fsusboot = boot.bootgenchunk_multimem(pi_fsus_sm, 200, 10, 1000, seed=3)
fsusboot = fsusboot.stack(year=("imem","isample"))
fsusboot['year'] = np.arange(0,2000,1)

amocboot = amocboot.isel(year=slice(0,pi_amoc.year.size))
fsusboot = fsusboot.isel(year=slice(0,pi_amoc.year.size))

# lagged regressions
lag = np.arange(-40,40)

min95_fsus = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='min95_fsusreg')
max95_fsus = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='max95_fsusreg')

bvals = xr.DataArray(np.zeros([len(lag), 1000]), coords=[lag, amocdif.iboot], 
                       dims=['lag','iboot'], name='b')

for ilag in np.arange(0,len(lag),1):
    amocshift = amocboot.shift(year=lag[ilag]).dropna('year')
    minyear = amocshift.year[0] ; maxyear = amocshift.year[amocshift.year.size-1]

    fsususe = fsusboot.where( (fsusboot.year >= minyear) & (fsusboot.year <= maxyear), drop=True)

    b = np.zeros([1000])
    for iboot in np.arange(0,1000,1):
        at, bt = linfit.linfit_xy(amocshift.isel(iboot=iboot), fsususe.isel(iboot=iboot))
        b[iboot] = bt
        bvals[ilag,iboot] = bt

    min95_fsus[ilag] = np.quantile(b, 0.025)
    max95_fsus[ilag] = np.quantile(b, 0.975)

# find the minimum regression coefficient and the associated lag
bmin = xr.DataArray(np.zeros([1000]), coords=[amocdif.iboot], dims=['iboot'], name='bmin')
lagmin = xr.DataArray(np.zeros([1000]), coords=[amocdif.iboot], dims=['iboot'], name='lagmin')

for iboot in np.arange(0,1000,1):
    bmin[iboot] = bvals.isel(iboot=iboot).min('lag')
    lagmin[iboot] = lag[np.argmin(np.array(bvals.isel(iboot=iboot)))]

#imin = np.argmin(bvals, 
#bvals = bvals.min(dim='lag')
#amocconstruct = bvals*amocdif
amocconstruct = xr.DataArray(np.zeros([1000, amocdif.year.size]), 
                   coords=[np.arange(0,1000,1),amocdif.year],
                   dims=['iboot','year'], name='amocconstruct') 
for iboot in np.arange(0,1000,1):
    laguse = lagmin.isel(iboot=iboot).values
    amocconstruct[iboot,:] = (bmin.isel(iboot=iboot).values)*\
                   amocdif.isel(iboot=iboot).shift(year=np.int(laguse))

#amocconstruct = bmin*amocdif.shift(year=lagmin

min95_bvals = bvals.quantile(0.025, dim='iboot')
max95_bvals = bvals.quantile(0.975, dim='iboot')

min95_amocconstruct = amocconstruct.quantile(0.025, dim='iboot')
max95_amocconstruct = amocconstruct.quantile(0.975, dim='iboot')

min95_amocconstruct = min95_amocconstruct.rename('min95_amocconstruct')
max95_amocconstruct = max95_amocconstruct.rename('max95_amocconstruct')


dat = xr.merge([min95_fsus, max95_fsus, min95_amocconstruct, max95_amocconstruct,
                min95_bvals, max95_bvals], compat='override')

dat.to_netcdf(pathout+'CESM2_FSUSregontoMOC.nc')



