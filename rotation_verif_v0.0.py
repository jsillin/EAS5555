'''
This script plots updraft_helicity from WRF over an hour of interest as well as
MRMS mid-level rotation for verification.

Written 3/29/21 by Jack Sillin for EAS 5555

This is version 0.0 of this script.
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
import scipy.ndimage as ndimage

#Read File
ncfile = Dataset("data/realtime_3-24-21_12z_333m/333mnest_v2_20z.nc")
mrmsfile = 'data/MRMS_03-25-20z_mlrotation.nc'

mrmsn = xr.open_dataset(mrmsfile).to_array().squeeze().isel(variable=0)
print(mrmsn)
# Get coordinates
lat = wrf.getvar(ncfile,'lat')
lon = wrf.getvar(ncfile,'lon')

lomax = np.max(lon).data
lomin = np.min(lon).data

mrms_lat_max = round(float(np.max(lat).values),2)
mrms_lat_min = round(float(np.min(lat).values),2)
mrms_lon_max = 360+round(float(lomax),2)
mrms_lon_min = 360+round(float(lomin),2)

lon_diff = mrms_lon_max-mrms_lon_min
lat_diff = mrms_lat_max-mrms_lat_min

num_lon_points = int(round(lon_diff/0.01,0))
num_lat_points = int(round(lat_diff/0.01,0))

lon_points = np.linspace(mrms_lon_min,mrms_lon_max,num_lon_points)
lat_points = np.linspace(mrms_lat_min,mrms_lat_max,num_lat_points)
mrms = mrmsn.isel(lon=(mrmsn.lon>mrms_lon_min)&(mrmsn.lon<mrms_lon_max))
mrms = mrms.reindex(lat=mrms.lat[::-1])
mrms = mrms.isel(lat=(mrms.lat>mrms_lat_min)&(mrms.lat<mrms_lat_max))
mrms_rotmax = mrms.max(dim='time')
print(mrms_rotmax)
print(np.max(mrms_rotmax))

uh_full_data = wrf.getvar(ncfile,'updraft_helicity',timeidx=wrf.ALL_TIMES)
print(uh_full_data)
uh_max_data = uh_full_data.max(dim='Time')
print(uh_max_data)
uh_max = ndimage.gaussian_filter(uh_max_data,sigma=3,order=0)
print(uh_max)

# Create the figure and axes
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111, projection = ccrs.PlateCarree())

# Add various geographical information to the plot
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)

uhmaxc = ax.contourf(lon,lat,uh_max,levels=range(50,1500,50),cmap='RdPu',extend='max')
mrmaxc = ax.contourf(lon_points,lat_points,mrms_rotmax, levels=range(5,40,5),cmap='Blues',extend='max')
time = wrf.getvar(ncfile,'times',timeidx=0)
dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

# Set titles
ax.set_title('Mid-Level Rotation Verification')
ax.set_title('Init 18z 3-23-21 from GFS',fontsize=11,loc='left')
ax.set_title('Valid: '+dtfs,fontsize=11,loc='right')
#ax1.set_extent((-88.1,-87.2,35,35.6))
ax.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')

# Save output graphic
plt.savefig('uhmax_verif_v2_333mnest_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
print('done')
'''
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
    wspd_10m = ((u_10m**2)+(v_10m**2))**.5

    #Make and print nice-looking datetime string
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    uh_filtered = ndimage.gaussian_filter(uh,sigma=3,order=0)
    srh_filtered = ndimage.gaussian_filter(srh,sigma=2,order=0)

    # Thin the wind arrays a bit for plotting nicely
    wind_slice = slice(20,-20,20)
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
    print(np.max(srh_filtered))
    # Plot 2m temps
    #capec = ax1.contourf(lon,lat,cape[0],levels=range(500,2500,100),alpha=0.5,extend='max',cmap='cool')
    srhc = ax1.contourf(lon,lat,srh_filtered,levels=range(100,2500,50),extend='both',cmap='cool')
    refc = ax1.contourf(lon,lat,simrad,norm=newnorm,cmap=newcmap,levels=range(5,75,1))
    uhc = ax1.contour(lon,lat,uh_filtered,colors='white',levels=range(100,1600,400))
    cbr = fig.colorbar(refc, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                        extendrect=False, ticks = range(-20,100,5), shrink=0.85)
    cbr.set_label('Reflectivity (dBZ)', fontsize = 14)

    #wsc = ax1.contour(lon,lat,wspd_10m,levels=range(10,50,5),cmap='PuRd',linewidths=1.5)
    # Plot 10m winds thinned out for better presentation
    ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

    # Set titles
    ax1.set_title('Simulated Reflectivity (dBZ), Storm-Relative Helicity (m^2/s^2) and 10m Winds (kts)')
    ax1.set_title('Init 18z 3-23-21 from GFS',fontsize=11,loc='left')
    ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')
    #ax1.set_extent((-88.1,-87.2,35,35.6))
    ax1.add_feature(USCOUNTIES.with_scale('5m'), edgecolor='dimgray')

    # Save output graphic
    plt.savefig('svrparams_v6_333mnest_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
'''
