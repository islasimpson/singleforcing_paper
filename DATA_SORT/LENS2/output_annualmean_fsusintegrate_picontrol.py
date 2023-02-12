import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read
from CASutils import modes_utils as modes

def calcglobalsum(darray, lon1, lon2, lat1, lat2):
    a = 6.37e6
    latrad = np.deg2rad(darray.lat)
    dphi = xr.DataArray(np.zeros([latrad.size]), coords=[darray.lat], dims=['lat'], name='dphi')
    dphimid = [ ((latrad.isel(lat=i+1).values - latrad.isel(lat=i).values)/2. +
                (latrad.isel(lat=i).values - latrad.isel(lat=i-1).values)/2.) for i in np.arange(1,latrad.size-1,1) ]
    dphi[1:darray.lat.size-1] = dphimid
    dphi[0] = (latrad.isel(lat=1).values - latrad.isel(lat=0).values)/2. + (latrad.isel(lat=0).values - np.deg2rad(-90))
    dphi[darray.lat.size-1] = (latrad.isel(lat=latrad.size-1).values - latrad.isel(lat=latrad.size-2).values)/2. + \
                               np.deg2rad(90) - latrad.isel(lat=latrad.size-1).values
    lonrad = np.deg2rad(darray.lon)
    dlon = lonrad[1].values - lonrad[0].values

    cosphi = xr.DataArray(np.zeros([latrad.size]), coords=[darray.lat], dims=['lat'], name='cosphi')
    cosphi[:] = np.cos(np.deg2rad(darray.lat.values))

    weights = xr.DataArray(np.ones([latrad.size, lonrad.size]), coords=[darray.lat, darray.lon], dims=['lat','lon'], name='weights')
    weights = a**2 * weights*cosphi*dphi*dlon

    region = darray.sel(lon=slice(lon1, lon2), lat=slice(lat1, lat2))
    weights = weights.sel(lon=slice(lon1, lon2), lat=slice(lat1,lat2))

    regionw = region.weighted(weights)
    integral = regionw.sum(("lon","lat"))

    return integral


topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/atm/month_1/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/"

fsnsam = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_FSNS_am.nc")

fsdsam = xr.open_dataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2/LENS2_FSDS_am.nc")

fsusam = fsdsam.FSDS - fsnsam.FSNS
fsusam = fsusam.rename('FSUS')

fsus_sum = calcglobalsum(fsusam, 0, 360, -90, 90)

fsus_sum.to_netcdf(pathout+"LENS2_fsus_integral_am.nc")



