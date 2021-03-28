'''
This script plots some severe weather parameters from a wrfout.nc file.
Intended to be used for storm-scale (several counties) visualizations

Written 3/28/21 by Jack Sillin for EAS 5555

This is version 0.0 of this script.

-Need to smooth and/or mask UH and SRH fields on the 333m nest data...
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
from metpy.plots import USCOUNTIES

#Read File
ncfile = Dataset("data/realtime_3-24-21_12z_333m/333mnest_v2_20z.nc")

wlon = -88.1
elon = -87.2
slat = 35
nlat = 35.6

#Loop though a lot of forecast times. It'll stop when you run out, don't worry.
for i in range(0,800):
    #Extract variables
    simrad = wrf.getvar(ncfile,'dbz',timeidx=i)
    slp = wrf.getvar(ncfile,'slp',timeidx=i)
    p = wrf.getvar(ncfile,'pressure',timeidx=i)
    #pwat = wrf.getvar(ncfile,'pw',timeidx=i)
    u = wrf.getvar(ncfile,'ua',timeidx=i)
    v = wrf.getvar(ncfile,'va',timeidx=i)
    u_10m,v_10m = wrf.getvar(ncfile,'uvmet10',timeidx=i)
    simrad=simrad.isel(bottom_top=2)
    time = wrf.getvar(ncfile,'times',timeidx=i)
    uh = wrf.getvar(ncfile,'updraft_helicity',timeidx=i)
    srh = wrf.getvar(ncfile, 'helicity',timeidx=i)
    vort = wrf.getvar(ncfile,'avo',timeidx=i)
    cape = wrf.getvar(ncfile,'cape_2d',timeidx=i)
    print(cape)
    wspd_10m = ((u_10m**2)+(v_10m**2))**.5
    print(np.max(wspd_10m))

    #Make and print nice-looking datetime string
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
    print(dtfs)
    print(uh)
    print(np.max(uh))
    print(srh)
    print(np.max(srh))
    print(vort)
    print(np.max(vort))

    # Get coordinates
    lat = wrf.getvar(ncfile,'lat')
    lon = wrf.getvar(ncfile,'lon')

    # Thin the wind arrays a bit for plotting nicely
    wind_slice = slice(9,-9,9)
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
    cbr = fig.colorbar(refc, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.85)
    cbr.set_label('Reflectivity (dBZ)', fontsize = 14)

    capec = ax1.contourf(lon,lat,cape[0],levels=range(500,2500,100),alpha=0.5,extend='max',cmap='cool')
    wsc = ax1.contour(lon,lat,wspd_10m,levels=range(10,50,5),cmap='PuRd',linewidths=1.5)
    # Plot 10m winds thinned out for better presentation
    ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

    # Set titles
    ax1.set_title('Simulated Reflectivity (dBZ), and 10m Winds (kts)')
    ax1.set_title('Init 18z 3-23-21 from GFS',fontsize=11,loc='left')
    ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')
    ax1.set_extent((-88.1,-87.2,35,35.6))
    ax1.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')

    # Save output graphic
    plt.savefig('svrparams_stormscale_v3_333mnest_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
