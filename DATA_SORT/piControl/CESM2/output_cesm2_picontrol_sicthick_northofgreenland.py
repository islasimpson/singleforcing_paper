import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import glob

from CASutils import averaging_utils as avg
from CASutils import lensread_utils as lens
from CASutils import calendar_utils as cal
from CASutils import readdata_utils as read

topdir="/project/mojave/cesm2/b.e21.B1850.f09_g17.CMIP6-piControl.001/ice/month_1/hi/"
pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/piControl/CESM2/"

files = glob.glob(topdir+"*.nc")

def fixtime(dat):
    timebnds = dat.time_bounds
    diff =  np.array(timebnds.isel(d2=1)) - np.array(timebnds.isel(d2=0))
    diff = diff/2.
    newtime = np.array(timebnds.isel(d2=0)) + diff
    dat['time'] = newtime
    return dat

dat = [ fixtime(xr.open_dataset(ifile, decode_times='False')).hi for ifile in files ]
dat = xr.concat(dat, dim='time')

dat = dat.groupby('time.year').mean('time')

tarea = xr.open_dataset(topdir+"b.e21.B1850.f09_g17.CMIP6-piControl.001.cice.h.hi.000101-009912.nc")
tarea = tarea.tarea

wgts = tarea
wgts = xr.where( (wgts.TLON > 240) & (wgts.TLON < 360), wgts, 0)
wgts = xr.where( (wgts.TLAT > 80) & (wgts.TLAT < 90), wgts, 0)

dat_w = dat.weighted(wgts.fillna(0))
dat_m = dat_w.mean(("ni","nj"))
dat_m.to_netcdf(pathout+"SICTHICK_annualmean_picontrol_northofgreenland.nc")
