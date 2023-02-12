import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from CASutils import averaging_utils as avg
from CASutils import bootstrap_utils as boot

landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig3/"


all1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1/"+\
                      "LENS1_TREFHT_am.nc")
all1['lon'] = landfrac.lon ; all1['lat'] = landfrac.lat
all1 = all1.sel(year=slice(1920,2050)).TREFHT

xghg1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
                        "XGHG_TREFHT_am.nc")
xghg1['lon'] = landfrac.lon ; xghg1['lat'] = landfrac.lat
xghg1 = xghg1.sel(year=slice(1920,2050)).TREFHT

xaer1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
                        "XAER_TREFHT_am.nc")
xaer1['lon'] = landfrac.lon ; xaer1['lat'] = landfrac.lat
xaer1 = xaer1.sel(year=slice(1920,2050)).TREFHT

xbmb1 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS1-SF/"+\
                        "XBMB_TREFHT_am.nc")
xbmb1['lon'] = landfrac.lon ; xbmb1['lat'] = landfrac.lat
xbmb1 = xbmb1.sel(year=slice(1920,2050)).TREFHT

all1gm = avg.cosweightlonlat(all1, 0, 360, -90, 90)
xghg1gm = avg.cosweightlonlat(xghg1, 0, 360, -90, 90)
xaer1gm = avg.cosweightlonlat(xaer1, 0, 360, -90, 90)
xbmb1gm = avg.cosweightlonlat(xbmb1, 0, 360, -90, 90)

bootall = boot.bootgen(all1gm, nsamples = all1.M.size, seed=3, nboots=1000)
bootxghg = boot.bootgen(xghg1gm, nsamples=xghg1.M.size, seed=4, nboots=1000)
bootxaer = boot.bootgen(xaer1gm, nsamples=xaer1.M.size, seed=5, nboots=1000)
bootxbmb = boot.bootgen(xbmb1gm, nsamples=xbmb1.M.size, seed=6, nboots=1000)

bootallm = bootall.mean('isample')
bootxghgm = bootxghg.mean('isample')
bootxaerm = bootxaer.mean('isample')
bootxbmbm = bootxbmb.mean('isample')

bootallmbase = bootallm.sel(year=slice(1920,1940)).mean('year')
bootxghgmbase = bootxghgm.sel(year=slice(1920,1940)).mean('year')
bootxaermbase = bootxaerm.sel(year=slice(1920,1940)).mean('year')
bootxbmbmbase = bootxbmbm.sel(year=slice(1920,1940)).mean('year')

bootallm = bootallm - bootallmbase
bootxghgm = bootxghgm - bootxghgmbase
bootxaerm = bootxaerm - bootxaermbase
bootxbmbm = bootxbmbm - bootxbmbmbase

bootghgm = bootallm - bootxghgm 
bootaerm = bootallm - bootxaerm
bootbmbm = bootallm - bootxbmbm
bootghgplusaerm = bootghgm + bootaerm
bootee = bootallm - (bootghgm + bootaerm + bootbmbm)

min95all = bootallm.quantile(0.025, dim='iboot') ; max95all = bootallm.quantile(0.975, dim='iboot')
min95ghg = bootghgm.quantile(0.025, dim='iboot') ; max95ghg = bootghgm.quantile(0.975, dim='iboot')
min95aer = bootaerm.quantile(0.025, dim='iboot') ; max95aer = bootaerm.quantile(0.975, dim='iboot')
min95bmb = bootbmbm.quantile(0.025, dim='iboot') ; max95bmb = bootbmbm.quantile(0.975, dim='iboot')
min95ee = bootee.quantile(0.025, dim='iboot') ; max95ee = bootee.quantile(0.975, dim='iboot')
min95ghgplusaer = bootghgplusaerm.quantile(0.025, dim='iboot')
max95ghgplusaer = bootghgplusaerm.quantile(0.975, dim='iboot')

min95all = min95all.rename('min95all') ; max95all = max95all.rename('max95all')
min95ghg = min95ghg.rename('min95ghg') ; max95ghg = max95ghg.rename('max95ghg')
min95aer = min95aer.rename('min95aer') ; max95aer = max95aer.rename('max95aer')
min95bmb = min95bmb.rename('min95bmb') ; max95bmb = max95bmb.rename('max95bmb')
min95ee = min95ee.rename('min95ee') ; max95ee = max95ee.rename('max95ee')
min95ghgplusaer = min95ghgplusaer.rename('min95ghgplusaer') 
max95ghgplusaer = max95ghgplusaer.rename('max95ghgplusaer')

all_ghg_aaer = xr.merge([min95all, max95all, min95ghg, max95ghg, min95aer, max95aer,
                         min95ghgplusaer, max95ghgplusaer], compat='override')
ee_bmb = xr.merge([min95bmb, max95bmb, min95ee, max95ee], compat='override') 

all_ghg_aaer.to_netcdf(pathout+'CESM1_TREFHTgm_boot_all_ghg_aaer.nc')
ee_bmb.to_netcdf(pathout+'CESM1_TREFHTgm_boot_ee_bmb.nc')
