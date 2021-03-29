'''
Verify WRF radar output against actual radar data

Written 3-25-21 by Jack Sillin for EAS 5555
Updated 3-28-21

This is version 0.1 of this script
-Added a colorbar (lol)
'''

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
import colorbars_v0 as cbars

mrmsn = xr.open_dataset('data/MRMS_BaseReflectivity_20210325_2200.grib2',engine='cfgrib')
#mrmsn = xr.open_dataset('https://thredds.ucar.edu/thredds/dodsC/grib/NCEP/MRMS/BaseRef/MRMS_BaseReflectivity_20210326_2100.grib2')
print(mrmsn)
mrms = mrmsn.isel(time=3).to_array().squeeze()
print(mrms)
time = mrms.valid_time
dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
print(dtfs)
mrmsn = mrms.where(mrms!=-999.9)
print(mrmsn)
print(np.max(mrmsn))
print(mrms)

#Read File
modelfile = Dataset("data/realtime_3-24-21_12z_1km/1kmnest.nc")

#Loop though a lot of forecast times. It'll stop when you run out, don't worry.
i=133
simrad = wrf.getvar(modelfile,'dbz',timeidx=i)
slp = wrf.getvar(modelfile,'slp',timeidx=i)
p = wrf.getvar(modelfile,'pressure',timeidx=i)
pwat = wrf.getvar(modelfile,'pw',timeidx=i)
u = wrf.getvar(modelfile,'ua',timeidx=i)
v = wrf.getvar(modelfile,'va',timeidx=i)
u_10m,v_10m = wrf.getvar(modelfile,'uvmet10',timeidx=i)
simrad=simrad.isel(bottom_top=2)
time = wrf.getvar(modelfile,'times',timeidx=i)

#Make and print nice-looking datetime string
dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())
print(dtfs)

# Get coordinates
lat = wrf.getvar(modelfile,'lat')
lon = wrf.getvar(modelfile,'lon')

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
mrms = mrmsn.isel(longitude=(mrmsn.longitude>mrms_lon_min)&(mrmsn.longitude<mrms_lon_max))
mrms = mrms.reindex(latitude=mrms.latitude[::-1])
mrms = mrms.isel(latitude=(mrms.latitude>mrms_lat_min)&(mrms.latitude<mrms_lat_max))
print(mrms)
print(np.max(mrms))
#mrms = mrms.filled(0)
#print(mrms)
#rint(np.max(mrms))
print(np.min(mrms))

# Thin the wind arrays a bit for plotting nicely
wind_slice = slice(15,-15,15)
sliced_lats = lat[wind_slice,wind_slice]
sliced_lons = lon[wind_slice,wind_slice]
sliced_u10m = u_10m[wind_slice,wind_slice]
sliced_v10m = v_10m[wind_slice,wind_slice]

# Create the figure and axes
fig = plt.figure(figsize=(20,12))
ax1 = fig.add_subplot(121, projection = ccrs.PlateCarree())

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

#h_contour = ax1.contour(lon, lat, slp, colors='dimgray', levels=range(940,1040,4),linewidths=2)

# Plot 10m winds thinned out for better presentation
#ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

# Set titles
ax1.set_title('1km WRF Simulated Reflectivity (dBZ)')
ax1.set_title('Init 18z 3-23-21 from GFS',fontsize=11,loc='left')
ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')



ax2 = fig.add_subplot(122, projection = ccrs.PlateCarree())

# Add various geographical information to the plot
ax2.coastlines(resolution='10m')
ax2.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
ax2.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
#ax2.set_extent((mrms_lon_min,mrms_lon_max,np.min(lat),np.max(lat)))

norm_ref, cmap_ref = ctables.registry.get_with_steps('NWSStormClearReflectivity', -20., .5)
newcmap = ListedColormap(cmap_ref(range(40, 194)))
newnorm = Normalize(0,77)
clevs = range(5,75,1)
# Plot 2m temps
refc = ax2.contourf(lon_points,lat_points,mrms.values,norm=newnorm,cmap=newcmap,levels=clevs)#,transform=ccrs.PlateCarree())
#tmp_2m32 = ax2.contour(lon,lat,t2m,colors='b', alpha = 0.8, levels = [32])


# Plot 10m winds thinned out for better presentation
#ax2.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

# Set titles
ax2.set_title('Observed Reflectivity (dBZ)')
ax2.set_title('MRMS',fontsize=11,loc='left')
ax2.set_title('Valid: '+dtfs,fontsize=11,loc='right')
fig.tight_layout()

cbars.addrefverifcolorbar(ax1,fig,refc,clevs)
fig.tight_layout()
plt.savefig('mrms_v20.png',bbox_inches='tight',pad_inches=0.1)
