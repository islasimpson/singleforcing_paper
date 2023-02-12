import xarray as xr
import numpy as np
import sys

avog = 6.022e23 # Avogadro's number
re = 6.3712e8 # Radius of the earth in cm

def convert_molecules_to_tg(dat,varname):
    """ Convert surface emissions in molecules/cm2/s to Tg 

    """
    # Start with moleculre/cm2/s.  Convert from molecules to grams.  
    # Divide by Avogadro's number to convert from molecules to moles.  
    #Multiply by molcular weight in g/mol to end up with g/cm2/s
    dat_g = dat[varname].molecular_weight*dat[varname]/avog

    # convert from per s to per year
    dat_g_y = dat_g*365.*86400.

    if "altitude" in dat.dims:
        print('you have altitudes')
        dz = dat.altitude_int[1:dat.altitude_int.size].values - dat.altitude_int[0:dat.altitude_int.size-1].values
        dz = xr.DataArray(dz, coords=[dat.altitude], dims=['altitude'], name='dz')
        dz = dz*1000.*100. # convert to cm
        dat_g_y = (dat_g_y*dz).sum('altitude')
    else if "lev" in dat[varname].dims:
        pressures = dat.hyam*1e5 + dat.hybm*dat.PS
        ipressures = dat.hyai*1e5 + dat.hybi*dat.PS
        dp = ipressures[1:ipressures.ilev.size,:,:,:] - np.array(ipressures[0:ipressures.ilev.size-1,:,:,:])
        dp = dp.rename({'ilev':'lev'})
        dp['lev'] = pressures.lev
        dz = (7.5/pressures)*dp # in km
        dat_g_y = (dat_g_y*dz).sum(lev)

    # Integrate over space
    dlon = np.deg2rad( (dat.lon[2] - dat.lon[1]))
    dlat = np.deg2rad( (dat.lat[2] - dat.lat[1]))
    area = xr.ones_like(dat.isel(time=0))
    weights = np.cos(np.deg2rad(area.lat))*dlat*dlon*re**2. # area in cm2
    dat_g_y_w = dat_g_y.weighted(weights)
    dattot = dat_g_y_w.sum(("lon","lat"))

    # Convert from grams to terra grams
    dattot = dattot/1e12

    return dattot

basepath="/project/cas/islas/python_savs/singleforcing_paper/EMISSIONS/CESM2/"

#----------------------------------------------------------------------------
#----------------------------------SO4---------------------------------------
#----------------------------------------------------------------------------

volc_so4_a1 = xr.open_dataset(basepath+'EE/SO4/ext/historical/'+\
'emissions-cmip6_so4_a1_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc')
volc_so4_a2 = xr.open_dataset(basepath+'EE/SO4/ext/historical/'+\
'emissions-cmip6_so4_a2_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc')

#----------------------------------------------------------------------------
#----------------------------------SO2---------------------------------------
#----------------------------------------------------------------------------
volc_so2 = xr.open_dataset(basepath+'EE/SO2/ext/historical/'+\
'emissions-cmip6_SO2_contvolcano_vertical_850-5000_0.9x1.25_c20170724.nc')



