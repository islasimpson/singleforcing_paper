import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys

from CASutils import bootstrap_utils as boot
from CASutils import averaging_utils as avg

landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/fig3/"

all2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"+\
   "LENS2_TREFHT_am.nc")
all2['lon'] = landfrac.lon ; all2['lat'] = landfrac.lat
all2 = all2.sel(year=slice(1920,2050)).TREFHT

ghg2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
   "GHG_TREFHT_am.nc")
ghg2['lon'] = landfrac.lon ; ghg2['lat'] = landfrac.lat
ghg2 = ghg2.sel(year=slice(1920,2050)).TREFHT

aer2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
   "AAER_TREFHT_am.nc")
aer2['lon'] = landfrac.lon ; aer2['lat'] = landfrac.lat
aer2 = aer2.sel(year=slice(1920,2050)).TREFHT

bmb2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
   "BMB_TREFHT_am.nc")
bmb2['lon'] = landfrac.lon ; bmb2['lat'] = landfrac.lat
bmb2 = bmb2.sel(year=slice(1920,2050)).TREFHT

ee2 = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
   "EE_TREFHT_am.nc")
ee2['lon'] = landfrac.lon ; ee2['lat'] = landfrac.lat
ee2 = ee2.sel(year=slice(1920,2050)).TREFHT

all2gm = avg.cosweightlonlat(all2, 0, 360, -90, 90)
ghg2gm = avg.cosweightlonlat(ghg2, 0, 360, -90, 90)
aer2gm = avg.cosweightlonlat(aer2, 0, 360, -90, 90)
bmb2gm = avg.cosweightlonlat(bmb2, 0, 360, -90, 90)
ee2gm = avg.cosweightlonlat(ee2, 0, 360, -90, 90)

bootall = boot.bootgen(all2gm, nsamples = all2gm.M.size, seed=3, nboots=1000)
bootghg = boot.bootgen(ghg2gm, nsamples = ghg2gm.M.size, seed=4, nboots=1000)
bootaer = boot.bootgen(aer2gm, nsamples = aer2gm.M.size, seed=5, nboots=1000)
bootbmb = boot.bootgen(bmb2gm, nsamples = bmb2gm.M.size, seed=6, nboots=1000)
bootee = boot.bootgen(ee2gm, nsamples = ee2gm.M.size, seed=7, nboots=1000)

bootallm = bootall.mean('isample')
bootghgm = bootghg.mean('isample')
bootaerm = bootaer.mean('isample')
bootbmbm = bootbmb.mean('isample')
booteem = bootee.mean('isample')

bootbaseall = bootallm.sel(year=slice(1920,1940)).mean('year')
bootbaseghg = bootghgm.sel(year=slice(1920,1940)).mean('year')
bootbaseaer = bootaerm.sel(year=slice(1920,1940)).mean('year')
bootbasebmb = bootbmbm.sel(year=slice(1920,1940)).mean('year')
bootbaseee = booteem.sel(year=slice(1920,1940)).mean('year')

bootallm = bootallm - bootbaseall
bootghgm = bootghgm - bootbaseghg
bootaerm = bootaerm - bootbaseaer
bootbmbm = bootbmbm - bootbasebmb
booteem = booteem - bootbaseee

bootghgplusaerm = bootghgm + bootaerm

min95all = bootallm.quantile(0.025, dim='iboot') ; max95all = bootallm.quantile(0.975, dim='iboot')
min95ghg = bootghgm.quantile(0.025, dim='iboot') ; max95ghg = bootghgm.quantile(0.975, dim='iboot')
min95aer = bootaerm.quantile(0.025, dim='iboot') ; max95aer = bootaerm.quantile(0.975, dim='iboot')
min95bmb = bootbmbm.quantile(0.025, dim='iboot') ; max95bmb = bootbmbm.quantile(0.975, dim='iboot')
min95ee = booteem.quantile(0.025, dim='iboot') ; max95ee = booteem.quantile(0.975, dim='iboot')
min95ghgplusaer = bootghgplusaerm.quantile(0.025, dim='iboot') 
max95ghgplusaer = bootghgplusaerm.quantile(0.975, dim='iboot')


min95all = min95all.rename('min95all') ; max95all = max95all.rename('max95all')
min95ghg = min95ghg.rename('min95ghg') ; max95ghg = max95ghg.rename('max95ghg')
min95aer = min95aer.rename('min95aer') ; max95aer = max95aer.rename('max95aer')
min95bmb = min95bmb.rename('min95bmb') ; max95bmb = max95bmb.rename('max95bmb')
min95ee = min95ee.rename('min95ee') ; max95ee = max95ee.rename('max95ee')
min95ghgplusaer = min95ghgplusaer.rename('min95ghgplusaer')
max95ghgplusaer = max95ghgplusaer.rename('max95ghgplusaer')

dat = xr.merge([min95all, max95all, min95ghg, max95ghg, min95aer, max95aer,
                min95bmb, max95bmb, min95ee, max95ee, min95ghgplusaer, max95ghgplusaer],
                compat='override')

dat.to_netcdf(pathout+'CESM2_TREFHTgm_boot.nc')






































