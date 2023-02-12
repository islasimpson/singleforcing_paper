import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from math import nan
from functools import partial

from CASutils import averaging_utils as avg
from CASutils import plotposition_utils as plotpos
from CASutils import colorbar_utils as cbars

import matplotlib as mpl
import sys

import warnings
warnings.filterwarnings('ignore')

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FOR_FIG8/"

plotpath="/project/cas/islas/python_plots/singleforcing_paper/figures_prelim/"
landfrac = xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

# pre-processor to ensure all lon and lat coordinates are the same.  Also adding flexibility for reading in seasonal data.
# Grabbing out 3 locations for analysis
lons=[286, 174, 95]
lats=[59, 68, 73]

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    try:
        year = ds.time.dt.year
        ds['time'] = year
        ds = ds.rename({'time':'year'})
    except:
        pass

    lons=[286, 174, 95]
    lats=[59, 68, 73]

    locations=[]
    for i in np.arange(0,len(lons),1):
        loc = ds.sel(lon=lons[i], lat=lats[i], method='nearest')
        locations.append(loc)
        
    locations = xr.concat(locations, dim='point')
    return locations

seas = 'jja' ; var = ['TREFHT','FSNO']

baselens2='/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/'
lens2 = []
for ivar in var:
    dat = xr.open_mfdataset(baselens2+'LENS2_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    lens2.append(dat)
lens2 = xr.merge(lens2).load()

basexaer2="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/CESM2-XAAER/"
xaer2 = []
for ivar in var:
    dat = xr.open_mfdataset(basexaer2+'xAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    xaer2.append(dat)
xaer2 = xr.merge(xaer2).load()

baseaer2="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"
aer2 = []
for ivar in var:
    dat = xr.open_mfdataset(baseaer2+'AAER_'+ivar+'_'+seas+'.nc', preprocess=partial(preprocessor))
    aer2.append(dat)
aer2 = xr.merge(aer2).load()

# Pick out the 1960-1980 anomalies

lens2_1960_1980 = lens2.sel(year=slice(1960,1980)).mean(['year','M'])
aer2_1960_1980 = aer2.sel(year=slice(1960,1980)).mean(['year','M'])
xaer2_1960_1980 = xaer2.sel(year=slice(1960,1980)).mean(['year','M'])

lens2_1920_1940 = lens2.sel(year=slice(1920,1940)).mean(['year','M'])
aer2_1920_1940 = aer2.sel(year=slice(1920,1940)).mean(['year','M'])
xaer2_1920_1940 = xaer2.sel(year=slice(1920,1940)).mean(['year','M'])

lens2anoms = lens2 - lens2_1920_1940
aer2anoms = aer2 - aer2_1920_1940
xaer2anoms = xaer2 - xaer2_1920_1940

lens2anoms_1960_1980 = lens2anoms.sel(year=slice(1960,1980)).mean(['year','M'])
aer2anoms_1960_1980 = aer2anoms.sel(year=slice(1960,1980)).mean(['year','M'])
xaer2anoms_1960_1980 = xaer2anoms.sel(year=slice(1960,1980)).mean(['year','M'])

lens2stack = lens2.stack(z=('year','M')).load()
aer2stack = aer2.stack(z=('year','M')).load()
xaer2stack = xaer2.stack(z=('year','M')).load()

lens2em = lens2.mean('M').load()
aer2em = aer2.mean('M').load()
xaer2em = xaer2.mean('M').load()


# Daily probability of no snow

fsnodaily = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FSNO/FSNO_JJA_locations_lens2_*.nc", 
                              combine='nested', concat_dim='M').FSNO
fsnodailyaer2 = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/FSNO/FSNO_JJA_locations_AAER_*.nc", 
                                  combine='nested', concat_dim='M').FSNO

years = lens2.year
pnosnow = np.zeros([years.size,len(lons)]) # probability of having snow

for iyear in np.arange(0,len(years),1):
    for iloc in np.arange(0,3,1):
        fsnouse = fsnodaily.isel(point=iloc)
        fsnouse = fsnouse.where( fsnouse.time.dt.year == years[iyear], drop=True)
        fsnouse = fsnouse.stack(z=("M","time"))
        
        ntime = fsnouse.z.size
        fnosno = fsnouse.where( fsnouse == 0, drop=True)
        pnosnow[iyear,iloc] = (fnosno.size/ntime)

years_aer2 = aer2.year
pnosnow_aer2 = np.zeros([years_aer2.size,len(lons)])
for iyear in np.arange(0,len(years_aer2),1):
    for iloc in np.arange(0,3,1):
        fsnouse = fsnodailyaer2.isel(point=iloc)
        fsnouse = fsnouse.where( fsnouse.time.dt.year == years[iyear] , drop=True)
        fsnouse = fsnouse.stack(z=("M","time"))
        
        ntime = fsnouse.z.size
        fnosno = fsnouse.where( fsnouse == 0, drop=True)
        pnosnow_aer2[iyear,iloc] = (fnosno.size/ntime)

pnosnow = xr.DataArray(pnosnow, coords=[years, np.arange(0,len(lons),1)],
              dims=['years', 'points'], name='pnosnow')

pnosnow_aer2 = xr.DataArray(pnosnow_aer2, coords=[years_aer2, np.arange(0,len(lons),1)],
                dims=['years','points'], name='pnosnow')



lens2stack = lens2stack.reset_index("z")
lens2em = lens2em

aer2stack = aer2stack.reset_index("z")
aer2em = aer2em

xaer2stack = xaer2stack.reset_index("z")
xaer2em = xaer2em


lens2stack.to_netcdf(pathout+'snow_lens2stack.nc')
lens2em.to_netcdf(pathout+'snow_lens2em.nc')
pnosnow.to_netcdf(pathout+'pnosnow_lens2.nc')

aer2stack.to_netcdf(pathout+'snow_aer2stack.nc')
aer2em.to_netcdf(pathout+'snow_aer2em.nc')
pnosnow_aer2.to_netcdf(pathout+'pnosnow_aer2.nc')

xaer2stack.to_netcdf(pathout+'snow_xaer2stack.nc')
xaer2em.to_netcdf(pathout+'snow_xaer2em.nc')







