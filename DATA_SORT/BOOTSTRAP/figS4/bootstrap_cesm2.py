import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
from functools import partial

from CASutils import bootstrap_utils as boot
from CASutils import linfit_utils as linfit
from CASutils import averaging_utils as avg

pathout="/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/BOOTSTRAP/figS4/"

varnames=['BURDENPOMdn', 'BURDENDUSTdn', 'BURDENSEASALTdn', 'BURDENSOAdn']

landfrac=xr.open_dataset("/project/cas/islas/cesm2le/fx/LANDFRAC_LENS2.nc")

def preprocessor(ds):
    ds['lon'] = landfrac.lon ; ds['lat'] = landfrac.lat
    ds = ds.sel(year=slice(1850,2050))
    return ds

for ivar in varnames:
    data = xr.open_mfdataset("/project/cas/islas/python_savs/singleforcing_paper/DATA_SORT/LENS2-SF/"+\
                "AAER_"+ivar+"_am.nc", preprocess=partial(preprocessor))

    datm = avg.cosweightlonlat(data[ivar], 0, 360, -90, 90)

    bootdat = boot.bootgen(datm, nsamples = data.M.size, seed=4, nboots=1000)
    bootdat3 = boot.bootgen(datm, nsamples = 3, seed=5, nboots=1000)

    bootdatm = bootdat.mean('isample')
    bootdat3m = bootdat3.mean('isample')

    min95 = bootdatm.quantile(0.025, dim='iboot') ; max95 = bootdatm.quantile(0.975, dim='iboot')
    min95_3 = bootdat3m.quantile(0.025, dim='iboot') ; max95_3 = bootdat3m.quantile(0.975, dim='iboot')

    min95 = min95.rename(ivar+'_min95') ; max95 = max95.rename(ivar+'_max95')
    min95_3 = min95_3.rename(ivar+'_min95_3') ; max95_3 = max95_3.rename(ivar+'_max95_3')

    if (ivar == varnames[0]):
        dat = xr.merge([min95, max95, min95_3, max95_3], compat='override')
    else:
        dat = xr.merge([dat, min95, max95, min95_3, max95_3], compat='override')

dat.to_netcdf(pathout+'CESM2_bootstrap.nc')





