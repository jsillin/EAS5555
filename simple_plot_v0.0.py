'''
This script makes a simple plot of surface temperature and 10m winds from a wrfout.nc file

Written 3/15/21 by Jack Sillin for EAS 5555

This is version 0.1 of this script.
'''

########## SETUP ##########

#Import Necessary Packages
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.lines as lines

# Read data with xarray
data_path = 'data/'
filename = 'matthew_wrfout.nc'
path = data_path + filename
print(path)
ds = xr.open_dataset(path)

########## PROCESSING DATA ##########

# Select/extract time of interest
ds = ds.isel(Time=1)
time = ds.XTIME
dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

# Get coordinates
lat = ds['XLAT']
lon = ds['XLONG']

# Extract 10m wind data and convert m/s to kts
u_10m = ds['U10']*1.94384449
v_10m = ds['V10']*1.94384449

# Extract 2m temp data and convert to F
t2m = ds['TH2']
t2m = ((t2m - 273.15)*(9./5.))+32.

# Thin the wind arrays a bit for plotting nicely
wind_slice = slice(3,-3,3)
sliced_lats = lat[wind_slice,wind_slice]
sliced_lons = lon[wind_slice,wind_slice]
sliced_u10m = u_10m[wind_slice,wind_slice]
sliced_v10m = v_10m[wind_slice,wind_slice]

########## PLOTTING ##########

# Define legend entries
blue_line = lines.Line2D([], [], color='b',label='32F Isotherm')

# Create the figure and axes
fig = plt.figure(figsize=(15,15))
ax1 = fig.add_subplot(111, projection = ccrs.PlateCarree())

# Add various geographical information to the plot
ax1.coastlines(resolution='10m')
ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
ax1.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)

# Plot 2m temps and 32F isotherm
tmp_2m = ax1.contourf(lon,lat,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(-20,100,5),transform=ccrs.PlateCarree())
tmp_2m32 = ax1.contour(lon,lat,t2m,colors='b', alpha = 0.8, levels = [32])
cbr = fig.colorbar(tmp_2m, orientation = 'horizontal', aspect = 80, ax = ax1, pad = 0.01,
                    extendrect=False, ticks = range(-20,100,5), shrink=0.7)
cbr.set_label('2m Temperature (F)', fontsize = 14)

# Plot 10m winds thinned out for better presentation
ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

# Set titles
ax1.set_title('2m Temperatures (F) and 10m Winds (kts)')
ax1.set_title('WRF Basic Case Study',fontsize=11,loc='right')
ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')

# Add legend
leg = ax1.legend(handles=[blue_line],loc=3,framealpha=1)

# Save output graphic
plt.savefig('sfc_tempsandwinds_v1'+dtfs+'.png')
