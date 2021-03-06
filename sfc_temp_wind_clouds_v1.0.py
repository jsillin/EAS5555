'''
This script makes a simple plot of surface temperature and 10m winds from a wrfout.nc file

Written 3/19/21 by Jack Sillin for EAS 5555
Updated 3/23/21

This is version 1.0 of this script.
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

# Read data with xarray
data_path = 'data/nested_matthew_v3/'
filename = 'matthew_innernest.nc'
path = data_path + filename
print(path)
dst = xr.open_dataset(path)

########## PROCESSING DATA ##########

# Select/extract time of interest

#print(ds['CLDFRA'])
#print(ds['RAINNC'])
#print(ds['PSFC']) #sfc pressure
#print(ds['P_HYD']) #hydrostatic pressure (MSLP?)
#print(ds['PHB']) #base-state geopotential
#print(ds['PH'])#perturbation geopotential
#print(ds['U'])#3d u wind
#print(ds['V'])#3d v wind
#print(ds['W'])#3d w wind

print('num_times:')
print(dst.XTIME)
for i in range(0,240):
    ds = dst.isel(Time=i)
    time = ds.XTIME
    print('num_times:')
    print(time)
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    # Get coordinates
    lat = ds['XLAT']
    lon = ds['XLONG']
    #print(lat.values)
    #print(lon.values)
    # Extract 10m wind data and convert m/s to kts
    u_10m = ds['U10']*1.94384449
    v_10m = ds['V10']*1.94384449

    # Extract 2m temp data and convert to F
    t2m = ds['TH2']
    t2m = ((t2m - 273.15)*(9./5.))+32.

    # Extract cloud fraction information
    cloud_frac = ds['CLDFRA']
    #print(np.shape(cloud_frac))
    total_cloud_frac = ds['CLDFRA'].max(dim='bottom_top')*100
    #print(total_cloud_frac)
    #print(t2m)
    #print(u_10m)
    #print(v_10m)
    '''
    new_ds = xr.Dataset(
    {
    "tot_cloud" : (["x","y","time"], total_cloud_frac),
    "2m_temp" : (["x","y","time"], t2m),
    "10m_uwind" : (["x","y","time"], u_10m),
    "10m_vwind" : (["x","y","time"], v_10m),
    },
    coords=
    {
    "lat":(["x","y"],lat.values),
    "lon":(["x","y"],lon.values),
    "time":time.values
    }
    )
    print(new_ds)
    '''
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
    tmp_2m = ax1.contourf(lon,lat,t2m,cmap='RdYlBu_r', alpha = 0.8, levels = range(60,90,2),transform=ccrs.PlateCarree())
    #tmp_2m32 = ax1.contour(lon,lat,t2m,colors='b', alpha = 0.8, levels = [32])
    cbars.addtempcolorbar(ax1,fig,tmp_2m,range(60,100,2))

    cldfrac = ax1.contourf(lon,lat,total_cloud_frac,cmap='bone',levels=range(0,101,1),alpha=0.3)
    cbars.addcloudcolorbar(ax1,fig,cldfrac,range(0,105,5))

    # Plot 10m winds thinned out for better presentation
    ax1.barbs(sliced_lons,sliced_lats,sliced_u10m,sliced_v10m, length=6,color='gray')

    # Set titles
    ax1.set_title('Total Cloud Cover 2m Temperatures (F) and 10m Winds (kts)')
    ax1.set_title('WRF Matthew Case Study',fontsize=11,loc='left')
    ax1.set_title('Valid: '+dtfs,fontsize=11,loc='right')

    # Add legend
    #leg = ax1.legend(handles=[blue_line],loc=3,framealpha=1)

    # Save output graphic
    plt.savefig('sfctempwind_cloudfrac_v4_innernest_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
