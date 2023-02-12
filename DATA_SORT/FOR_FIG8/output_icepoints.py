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

# Grabbing out 3 locations for analysis
lons=[215, 85, 275]
lats=[71, 79, 60]

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    try:
        year = ds.time.dt.year
        ds['time'] = year
        ds = ds.rename({'time':'year'})
    except:
        pass

    lons=[215, 85, 275]
    lats=[71, 79, 60]

    locations=[]
    for i in np.arange(0,len(lons),1):
        loc = ds.sel(lon=lons[i], lat=lats[i], method='nearest')
        locations.append(loc)

    locations = xr.concat(locations, dim='point')
    return locations

seas = 'jja' ; var = ['TREFHT','ICEFRAC']

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

# Daily probability of zero sea ice

icefracdaily = xr.open_mfdataset(
             "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/aice_d/*lens2*.nc",
               combine='nested', concat_dim='M').aice_d
icefracdailyaer2 = xr.open_mfdataset(
             "/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/aice_d/*aaer2*.nc",
               combine='nested', concat_dim='M').aice_d


years = lens2.year
pnoice = np.zeros([years.size, len(lons)]) # probability of having no ice
for iyear in np.arange(0,len(years),1):
    for iloc in np.arange(0,3,1):
        aiceuse = icefracdaily.isel(point=iloc)
        aiceuse = aiceuse.where( aiceuse.time.dt.year == years[iyear], drop=True)
        aiceuse = aiceuse.stack(z=("M","time"))
   
        ntime = aiceuse.z.size
        noice = aiceuse.where( aiceuse == 0, drop=True)
        pnoice[iyear,iloc] = (noice.size/ntime)

years_aer2 = aer2.year
pnoice_aer2 = np.zeros([years_aer2.size, len(lons)])
for iyear in np.arange(0,len(years_aer2),1):
    for iloc in np.arange(0,3,1):
        aiceuse = icefracdailyaer2.isel(point=iloc)
        aiceuse = aiceuse.where( aiceuse.time.dt.year == years[iyear], drop=True)
        aiceuse = aiceuse.stack(z=("M","time"))

        ntime = aiceuse.z.size
        noice = aiceuse.where( aiceuse == 0, drop=True)
        pnoice_aer2[iyear, iloc] = (noice.size/ntime)

pnoice = xr.DataArray(pnoice, coords=[years, np.arange(0,len(lons),1)],
                      dims=['years','points'], name='pnoice')
pnoice_aer2 = xr.DataArray(pnoice_aer2, coords=[years_aer2, np.arange(0,len(lons),1)],
                      dims=['years','points'], name='pnoice')

lens2stack = lens2stack.reset_index("z")
lens2em = lens2em

aer2stack = aer2stack.reset_index("z")
aer2em = aer2em

xaer2stack = xaer2stack.reset_index("z")
xaer2em = xaer2em

lens2stack.to_netcdf(pathout+'ice_lens2stack.nc')
lens2em.to_netcdf(pathout+'ice_lens2em.nc')
pnoice.to_netcdf(pathout+'pnoice_lens2.nc')

aer2stack.to_netcdf(pathout+'ice_aer2stack.nc')
aer2em.to_netcdf(pathout+'ice_aer2em.nc')
pnoice_aer2.to_netcdf(pathout+'pnoice_aer2.nc')

xaer2stack.to_netcdf(pathout+'ice_xaer2stack.nc')
xaer2em.to_netcdf(pathout+'ice_xaer2em.nc')
    























