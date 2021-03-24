'''
This script makes a plots MSLP, simulated reflectivity, and 10m winds from a wrfout.nc file

Written 3/19/21 by Jack Sillin for EAS 5555
Updated 3/24/21

This is version 0.2 of this script.
'''

########## SETUP ##########

#Import Necessary Packages
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.lines as lines
import colorbars_v0 as cbars
import wrf
from netCDF4 import Dataset
from metpy.plots import ctables
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize

#Read File
ncfile = Dataset("data/realtime_3-23-21_18z_1km/middlenest3km.nc")

#Loop though a lot of forecast times. It'll stop when you run out, don't worry.
for i in range(0,240):
    #Extract variables
    simrad = wrf.getvar(ncfile,'dbz',timeidx=i)
    slp = wrf.getvar(ncfile,'slp',timeidx=i)
    p = wrf.getvar(ncfile,'pressure',timeidx=i)
    pwat = wrf.getvar(ncfile,'pw',timeidx=i)
    u = wrf.getvar(ncfile,'ua',timeidx=i)
    v = wrf.getvar(ncfile,'va',timeidx=i)
    u_10m,v_10m = wrf.getvar(ncfile,'uvmet10',timeidx=i)
    simrad=simrad.isel(bottom_top=2)
    time = wrf.getvar(ncfile,'times',timeidx=i)

    #Make and print nice-looking datetime string
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    print(dtfs)

    # Get coordinates
    lat = wrf.getvar(ncfile,'lat')
    lon = wrf.getvar(ncfile,'lon')

    # Thin the wind arrays a bit for plotting nicely
    wind_slice = slice(6,-6,6)
    sliced_lats = lat[wind_slice,wind_slice]
    sliced_lons = lon[wind_slice,wind_slice]
    sliced_u10m = u_10m[wind_slice,wind_slice]
    sliced_v10m = v_10m[wind_slice,wind_slice]

    # Create the figure and axes
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = ccrs.PlateCarree())

    # Add various geographical information to the plot
    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
    ax1.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)

    norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
    newcmap = ListedColormap(cmap_ref(range(40, 194)))
    newnorm = Normalize(0,77)

    # Plot 2m temps
    refc = ax1.contourf(lon,lat,simrad,norm=newnorm,cmap=newcmap,levels=range(5,75,1))
    #tmp_2m32 = ax1.contour(lon,lat,t2m,colors='b', alpha = 0.8, levels = [32])
    cbr = fig.colorbar(refc, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.85)
    cbr.set_label('Reflectivity (dBZ)', fontsize = 14)

    h_contour = ax1.contour(lon, lat, slp, colors='dimgray', levels=range(940,1040,4),linewidths=2)

    # Plot 10m winds thinned out for better presentation
    ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

    # Set titles
    ax1.set_title('Simulated Reflectivity (dBZ), MSLP (hPa), and 10m Winds (kts)')
    ax1.set_title('Init 18z 3-23-21 from GFS',fontsize=11,loc='left')
    ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')

    # Add legend
    #leg = ax1.legend(handles=[blue_line],loc=3,framealpha=1)

    # Save output graphic
    plt.savefig('simrad_realtime_v1_3kmnest_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
