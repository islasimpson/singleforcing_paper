import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

from CASutils import bootstrap_utils as boot

#---------------------------------------------------------------------------
#--Bootstrapping land and not land averaged surface upward SW for uncertainty
#--range on Figure 6
#--Data pre-processed at 
#--~/DATA_SORT/XXXX/output_allmonths_50n90n.py
#--where XXXX is the experiment
#---------------------------------------------------------------------------

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig6/"

baseaer2='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/'
fsds_aer2 = xr.open_dataset(baseaer2+'AAER_FSDS_50n90n_land_notland.nc')
fsns_aer2 = xr.open_dataset(baseaer2+'AAER_FSNS_50n90n_land_notland.nc')
fsus_aer2 = fsds_aer2 - fsns_aer2

year = fsus_aer2.land.isel(M=0)[fsus_aer2.time.dt.month == 1].time.dt.year
fsus_aer2_monthly_land = xr.DataArray(np.zeros([fsus_aer2.M.size, int(fsus_aer2.time.size/12), 12]),
                                        dims=['M','year','month'], 
                                        coords=[fsus_aer2.M, year, np.arange(0,12,1)], name='FSUS')
fsus_aer2_monthly_notland = xr.DataArray(np.zeros([fsus_aer2.M.size, int(fsus_aer2.time.size/12), 12]),
                                        dims=['M','year','month'], 
                                        coords=[fsus_aer2.M, year, np.arange(0,12,1)], name='FSUS')

for imon in np.arange(1,12,1):
    fsus_aer2_monthly_land[:,:,imon-1] = fsus_aer2.land[:,fsus_aer2.time.dt.month == imon]
    fsus_aer2_monthly_notland[:,:,imon-1] = fsus_aer2.notland[:,fsus_aer2.time.dt.month == imon]

fsus_aer2_land_base = fsus_aer2_monthly_land.sel(year=slice(1920,1940)).mean('year')
fsus_aer2_notland_base = fsus_aer2_monthly_notland.sel(year=slice(1920,1940)).mean('year')

fsus_aer2_monthly_land = fsus_aer2_monthly_land - fsus_aer2_land_base
fsus_aer2_monthly_notland = fsus_aer2_monthly_notland - fsus_aer2_notland_base

def calc21ymeans(dat):
    datm = dat.rolling(year=21, min_periods=21, center='True').mean('year').dropna('year')
    return datm

fsus_aer2_monthly_land_21y = calc21ymeans(fsus_aer2_monthly_land)
fsus_aer2_monthly_notland_21y = calc21ymeans(fsus_aer2_monthly_notland)

boot_land = boot.bootgen(fsus_aer2_monthly_land_21y, nsamples=3, seed=5, nboots=1000)
boot_notland = boot.bootgen(fsus_aer2_monthly_notland_21y, nsamples=3, seed=6, nboots=1000)

boot_land = boot_land.mean('isample')
boot_notland = boot_notland.mean('isample')

boot_land = boot_land - fsus_aer2_monthly_land.mean('M')
boot_notland = boot_notland - fsus_aer2_monthly_notland.mean('M')

min95_land = boot_land.quantile(0.025, dim='iboot') 
max95_land = boot_land.quantile(0.975, dim='iboot')
min95_notland = boot_notland.quantile(0.025, dim='iboot') 
max95_notland = boot_notland.quantile(0.975, dim='iboot')

min95_land = min95_land.rename('min95_land')
max95_land = max95_land.rename('max95_land')
min95_notland = min95_notland.rename('min95_notland')
max95_notland = max95_notland.rename('max95_notland')

dat = xr.merge([min95_land, max95_land, min95_notland, max95_notland], compat='override')

dat.to_netcdf(pathout+'CESM2_bootstrap.nc')

