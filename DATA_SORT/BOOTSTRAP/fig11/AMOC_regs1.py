
import sys
from CASutils import linfit_utils as linfit
from CASutils import bootstrap_utils as boot

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig11/"

pi_amoc1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing/prelim_jun2022/danabasogluetal2019/outputdata/CESM1/piControl_AMOC_am.nc")
pi_amoc1 = pi_amoc1.MOC

pi_rho1 = xr.open_dataset(
    "/project/cas/islas/python_savs/singleforcing/prelim_jun2022/danabasogluetal2019/outputdata/CESM1/cesm1_picontrol_rho_top203_60to35W_50to65N.nc")
pi_rho1 = pi_rho1.RHO

pi_rhoparts1 = xr.open_dataset(
    "/project/cas/islas/python_savs/singleforcing/prelim_jun2022/danabasogluetal2019/outputdata/CESM1/RHO_contributions_CESM1_picontrol.nc")

pi_temp1 = pi_rhoparts1.TEMP_part[pi_rhoparts1.time.dt.month == 3]
pi_salt1 = pi_rhoparts1.SALT_part[pi_rhoparts1.time.dt.month == 3]
pi_rhofromeos1 = pi_rhoparts1.RHO_from_EOS[pi_rhoparts1.time.dt.month == 3]

pi_amoc1anoms = pi_amoc1 - pi_amoc1.mean('year')
pi_rho1anoms = pi_rho1 - pi_rho1.mean('time')
pi_rhofromeos1anoms = pi_rhofromeos1 - pi_rhofromeos1.mean('time')
pi_temp1anoms = pi_temp1 - pi_temp1.mean('time')
pi_salt1anoms = pi_salt1 - pi_salt1.mean('time')


# 10 year running means
pi_amoc1_sm = pi_amoc1anoms.rolling(year=10, min_periods=10, center='True').mean('year').dropna('year')
pi_rhofromeos1_sm = pi_rhofromeos1anoms.rolling(time=10, min_periods=10, center='True').mean('time').dropna('time')
pi_salt1_sm = pi_salt1anoms.rolling(time=10, min_periods=10, center='True').mean('time').dropna('time')
pi_temp1_sm = pi_temp1anoms.rolling(time=10, min_periods=10, center='True').mean('time').dropna('time')


amocboot = boot.bootgenchunk_multimem(pi_amoc1_sm, 200, 9, 1000, seed=3)
amocboot = amocboot.stack(year=("imem","isample"))
amocboot['year'] = np.arange(0,1800,1)

rhofromeosboot = boot.bootgenchunk_multimem(pi_rhofromeos1_sm, 200, 9, 1000, seed=3)
rhofromeosboot = rhofromeosboot.stack(year=("imem","isample"))
rhofromeosboot['year'] = np.arange(0,1800,1)

saltboot = boot.bootgenchunk_multimem(pi_salt1_sm, 200, 9, 1000, seed=3)
saltboot = saltboot.stack(year=("imem","isample"))
saltboot['year'] = np.arange(0,1800,1)

tempboot = boot.bootgenchunk_multimem(pi_temp1_sm, 200, 9, 1000, seed=3)
tempboot = tempboot.stack(year=("imem","isample"))
tempboot['year'] = np.arange(0,1800,1)

amocboot = amocboot.isel(year=slice(0,pi_amoc1_sm.year.size))
rhofromeosboot = rhofromeosboot.isel(year=slice(0,pi_amoc1_sm.year.size))
saltboot = saltboot.isel(year=slice(0,pi_amoc1_sm.year.size))
tempboot = tempboot.isel(year=slice(0,pi_amoc1_sm.year.size))



# lagged regressions
lag = np.arange(-40,40)

min95_amoc = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='min95_amoc')
max95_amoc = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='max95_amoc')

min95_rho = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='min95_rho')
max95_rho = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='max95_rho')

min95_salt = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='min95_salt')
max95_salt = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='max95_salt')

min95_temp = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='min95_temp')
max95_temp = xr.DataArray(np.zeros([len(lag)]), coords=[lag], dims=['lag'], name='max95_temp')


for ilag in np.arange(0,len(lag),1):

    #---------CESM1
    amocshift = amocboot.shift(year=lag[ilag]).dropna('year')
    minyear = amocshift.year[0] ; maxyear = amocshift.year[amocshift.year.size-1]

    amocuse = amocboot.where( (amocboot.year >= minyear) & (amocboot.year <= maxyear), drop=True ) 
  
    b = np.zeros([1000])
    for iboot in np.arange(0,1000,1):
        at, bt = linfit.linfit_xy(amocshift.isel(iboot=iboot), amocuse.isel(iboot=iboot))
        b[iboot] = bt

    min95_amoc[ilag] = np.quantile(b, 0.025)
    max95_amoc[ilag] = np.quantile(b, 0.975)

    rhouse = rhofromeosboot.where( (rhofromeosboot.year >= minyear) & (rhofromeosboot.year <= maxyear), 
                  drop=True)
    b = np.zeros([1000])
    for iboot in np.arange(0,1000,1):
        at, bt = linfit.linfit_xy(amocshift.isel(iboot=iboot), rhouse.isel(iboot=iboot))
        b[iboot] = bt

    min95_rho[ilag] = np.quantile(b, 0.025)
    max95_rho[ilag] = np.quantile(b, 0.975)

    saltuse = saltboot.where( (saltboot.year >= minyear) & (saltboot.year <= maxyear), 
                  drop=True)
    b = np.zeros([1000])
    for iboot in np.arange(0,1000,1):
        at, bt = linfit.linfit_xy(amocshift.isel(iboot=iboot), saltuse.isel(iboot=iboot))
        b[iboot] = bt

    min95_salt[ilag] = np.quantile(b, 0.025)
    max95_salt[ilag] = np.quantile(b, 0.975)

    tempuse = tempboot.where( (tempboot.year >= minyear) & (tempboot.year <= maxyear), 
                  drop=True)
    b = np.zeros([1000])
    for iboot in np.arange(0,1000,1):
        at, bt = linfit.linfit_xy(amocshift.isel(iboot=iboot), tempuse.isel(iboot=iboot))
        b[iboot] = bt

    min95_temp[ilag] = np.quantile(b, 0.025)
    max95_temp[ilag] = np.quantile(b, 0.975)


dat = xr.merge([min95_amoc, max95_amoc,
                min95_rho, max95_rho,
                min95_salt, max95_salt,
                min95_temp, max95_temp], compat='override')

dat.to_netcdf(pathout+'CESM1_MOCregressions.nc')
