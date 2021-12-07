# An interface for auto-plotting basic soundingmaps
# Written by Jack Sillin, last updated 3/6/21 version 0.2

# Imports
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import xarray as xr
import metpy
from datetime import datetime
import datetime as dt
from metpy.units import units
import scipy.ndimage as ndimage
from metpy.plots import USCOUNTIES
import cartopy
from scipy.ndimage.filters import generic_filter as gf
from metpy.plots import USCOUNTIES
from metpy.plots import SkewT
import metpy.calc as mpcalc
from math import log, exp
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import supplementary_tools as spt
import soundingmaps as smap

#Set interactive mode on (True) or off (False)
interactive = True

#In interactive mode, user is prompted for domain info
if interactive == True:
    # Welcome message for user
    print("\nWelcome to Jack Sillin's Soundingmap Plotting Tool!")
    print('\n \nNow you can make your own sounding maps with just a couple inputs:')
    print('\nDomain Info')
    print("Domain Size: string either 'local' (wfo-size) or 'regional' (a bit bigger)")
    print('Domain Name: whatever you want (string)')
    print('Center Lat: latitude of the centerpoint of your domain (float) in deg N')
    print('Center Lon: longitude of the centerpoint of your domain (flotat) in deg W')
    print('\nSounding Plot Info:')
    print("Want to shade CAPE/CIN on soundings? Options are 'yes' or 'no'")
    print("Want to plot wet bulb profiles? Options are 'yes' or 'no'")
    print('\nOther: ')
    print("Model: forecast model to plot (string). Options are 'GFS','NAM' ('RAP' at your own risk)")
    print("What phenomenon are you interested in? Options are 'winter' or 'severe'")
    print('\n \nMore features and options will be added soon. Happy plotting!')

    # User-Defined Inputs
    domainsize=input('Domain Size: ')
    domainname=input('Domain Name: ')
    centerlat =float(input('Center Lat: '))
    centerlon =float(input('Center Lon: '))

    cape_input = input('Shade CAPE? ')
    wetbulb_input = input('Plot Wet Bulb? ')

    model = input('Model: ')
    season = input('Season: ')

# In quiet mode, domain info is set with variables here
else:
    domainsize='regional'
    domainname='raleigh'
    centerlat = 35.8
    centerlon = 78.8
    model = 'GFS'
    cape_input = 'no'
    wetbulb_input = 'yes'
    season = 'winter'

# Regardless of which mode is used, a few more variables will need to be set
# based on the domain of interest:
if domainsize=='regional':
    south = centerlat-6
    north = centerlat+6
    east = 360-(centerlon-7.5)
    west = 360-(centerlon+7.5)

elif domainsize=='local':
    south = centerlat-1.625
    north = centerlat+1.625
    east = 360-(centerlon-2)
    west = 360-(centerlon+2)

# Convert CAPE/Wetbulb inputs to booleans
if cape_input == 'yes':
    cape_bool = True
elif cape_input == 'no':
    cape_bool = False
else:
    print('Invalid Input For CAPE (Need to pick "yes" or "no")')

if wetbulb_input == 'yes':
    wetbulb_bool = True
elif wetbulb_input == 'no':
    wetbulb_bool = False
else:
    print('Invalid Input For wetbulb (Need to pick "yes" or "no")')

#Prepare to pull data from NOMADS
url = spt.get_url(model)
mdate, init_hour = spt.get_init_time(model)

# Create new directory
output_dir = mdate+'_'+init_hour+'00'
spt.mkdir_p(output_dir)
spt.mkdir_p(output_dir+'/'+model)

#Parse data using MetPy
ds = xr.open_dataset(url)
times = ds['tmp2m'].metpy.time
init_time = ds['time'][0]

# If plotting GFS data, subset to North America to save some time
if model == 'GFS':
    lats = np.arange(20,55,0.25)
    lons = np.arange(240,310,0.25)
    ds = ds.sel(lat=lats,lon=lons)

#Get timestep info based on which model you're plotting.
timestepinfo = spt.get_num_timesteps(model)
etime = timestepinfo[0]
delt = timestepinfo[1]

if season == 'winter':
    #Run through each forecast hour and make the plots
    for i in range(0,etime):
        #Get the data for the forecast hour of interest
        data = ds.metpy.parse_cf()
        data = data.isel(time=i)

        #Rename variables to useful things
        data = data.rename(spt.get_varlist(model))

        #Pull out the categorical precip type arrays
        catrain = data['catrain'].squeeze()
        catsnow = data['catsnow'].squeeze()
        catsleet = data['catsleet'].squeeze()
        catice = data['catice'].squeeze()

        cape = data['cape'].squeeze()
        u10m = data['u'].squeeze()
        v10m = data['v'].squeeze()
        u10m = u10m*1.94384449
        v10m = v10m*1.94384449
        wspd = ((u10m**2)+(v10m**2))**.5

        #This extends each ptype one gridpoint outwards to prevent a gap between
        #different ptypes
        radius = 1
        kernel = np.zeros((2*radius+1,2*radius+1))
        y1,x1 = np.ogrid[-radius:radius+1,-radius:radius+1]
        mask=x1**2+y1**2 <=radius**2
        kernel[mask]=1

        #Make the ptype arrays nicer looking
        snowc= gf(catsnow,np.max,footprint=kernel)
        icec = gf(catice,np.max,footprint=kernel)
        sleetc = gf(catsleet,np.max,footprint=kernel)
        rainc = gf(catrain,np.max,footprint=kernel)

        #Coordinate stuff
        #vertical, = data['temperature'].metpy.coordinates('vertical')
        time = data['temperature'].metpy.time
        x, y = data['temperature'].metpy.coordinates('x', 'y')
        lats, lons = xr.broadcast(y, x)
        zH5_crs = data['temperature'].metpy.cartopy_crs
        wind_slice = spt.get_windslice(model,domainsize)
        #Processing surface temperature data
        t2m = data['sfc_temp'].squeeze()
        t2m = ((t2m - 273.15)*(9./5.))+32.

        td2m = data['sfc_td'].squeeze()
        td2m = ((td2m - 273.15)*(9./5.))+32.
        td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

        #Fetch reflectivity data
        reflectivity = data['radar'].squeeze()

        #Create masked arrays for each ptype
        rain = np.ma.masked_where(rainc==0,reflectivity)
        sleet = np.ma.masked_where(sleetc==0,reflectivity)
        ice = np.ma.masked_where(icec==0,reflectivity)
        snow = np.ma.masked_where(snowc==0,reflectivity)

        #Process pressure level temp/rh data
        prs_temps = data['temperature']
        prs_relh = data['rh']

        #Process MSLP data
        mslp = data['mslp']/100.
        mslpc = mslp.squeeze()
        mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)
        sfc_pressure = data['spres'].squeeze()/100

        #This creates a nice-looking datetime label
        dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

        ########## SET UP FIGURE ##################################################
        fig = plt.figure(figsize=(15,15))
        ax1 = fig.add_subplot(111, projection = zH5_crs)

        ax1.coastlines(resolution='10m')
        ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
        ax1.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
        ax1.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='silver')

        ########## PLOTTING #######################################################
        #Plot 2m 32F isotherm
        tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])

        #Plot labeled MSLP contours
        h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1080,4),linewidths=1,alpha=0.7)
        #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

        #Define levels and colors for plotting precip types
        ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
        qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
        qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
        qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
        qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']

        #Plot surface-based CAPE
        capep = ax1.contourf(x, y, cape, levels=[100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000], alpha = 0.6, cmap='RdPu')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])

        #Plot the underlying precip type shadings
        #Use try/except so that if no ice/sleet/etc is present, things don't break
        try:
            ra = ax1.contourf(x,y,rain,colors=qr_cols,levels=ref_levs,alpha=0.4,extend='max')
        except:
            print('no rain')
        try:
            sn = ax1.contourf(x,y,snow,colors=qs_cols,levels=ref_levs,alpha=0.4,extend='max')
        except:
            print('no snow')
        try:
            ip = ax1.contourf(x,y,sleet,colors=qi_cols,levels=ref_levs,alpha=0.4,extend='max')
        except:
            print('no sleet')
        try:
            zr = ax1.contourf(x,y,ice, colors=qz_cols,levels=ref_levs,alpha=0.4,extend='max')
        except:
            print('no ice')

        #Plot 10m wind barbs
        ax1.barbs(x[wind_slice],y[wind_slice],u10m[wind_slice,wind_slice],v10m[wind_slice,wind_slice], length=6,color='gray')

        #Plot soundings
        smap.plot_soundings(fig,ax1,prs_temps,prs_relh,sfc_pressure,centerlat,centerlon,domainsize,model,cape=cape_bool,wetbulb=wetbulb_bool)
        #Set plot extent and titles
        ax1.set_extent((west,east,south,north))
        ax1.set_title('Precipitation Type, MSLP, and Selected Soundings',fontsize=16)
        ax1.set_title('Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
        ax1.set_title(model+' Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

        #Save plot
        plt.savefig(output_dir+'/'+model+'/'+domainname+'_winter_soundingmap'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)

elif season=='severe':
    #Run through each forecast hour and make the plots
    for i in range(0,etime):
        #Get the data for the forecast hour of interest

        #Parse to get crs via metpy:

        #This is the old way for metpy 0.12.2
        #data = ds.metpy.parse_cf()

        #This is the new way for metpy 1.0.0
        data = ds.metpy.assign_crs(grid_mapping_name='latitude_longitude')

        #Select forecast hour of interest
        data = data.isel(time=i)

        #Rename variables to useful things
        data = data.rename(spt.get_varlist(model))

        #Extract some data into variable-specific arrays
        cape = data['cape'].squeeze()
        u10m = data['u'].squeeze()
        v10m = data['v'].squeeze()
        u10m = u10m*1.94384449
        v10m = v10m*1.94384449
        wspd = ((u10m**2)+(v10m**2))**.5
        sfc_pressure = data['spres'].squeeze()/100

        #Coordinate stuff
        #vertical, = data['temperature'].metpy.coordinates('vertical')
        time = data['temperature'].metpy.time
        x, y = data['temperature'].metpy.coordinates('x', 'y')
        lats, lons = xr.broadcast(y, x)
        zH5_crs = data['temperature'].metpy.cartopy_crs
        wind_slice = spt.get_windslice(model,domainsize)

        #Processing surface temperature/dew point data
        t2m = data['sfc_temp'].squeeze()
        t2m = ((t2m - 273.15)*(9./5.))+32.
        td2m = data['sfc_td'].squeeze()
        td2m = ((td2m - 273.15)*(9./5.))+32.
        td2ms = ndimage.gaussian_filter(td2m,sigma=5,order=0)

        #Fetch reflectivity data
        reflectivity = data['radar'].squeeze()

        #Process pressure level temp/rh data
        if model =='GFS':
            ilev = 18
        elif model == 'NAM':
            ilev = 30
        elif model == 'RAP':
            ilev = 24

        prs_temps = data['temperature'].isel(lev=slice(1,ilev,1))
        prs_relh = data['rh'].isel(lev=slice(1,ilev,1))
        print(prs_temps)
        #Process MSLP data
        mslp = data['mslp']/100.
        mslpc = mslp.squeeze()
        mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)

        #This creates a nice-looking datetime label
        dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

        ########## SET UP FIGURE ##################################################
        fig = plt.figure(figsize=(15,15))
        ax1 = fig.add_subplot(111, projection = zH5_crs)

        #Add various geographical information to the plot
        ax1.coastlines(resolution='10m')
        ax1.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=1.5)
        ax1.add_feature(cfeature.STATES.with_scale('10m'), linewidth=2.0)
        ax1.add_feature(USCOUNTIES.with_scale('500k'), edgecolor='silver')

        ########## PLOTTING #######################################################
        #Plot 2m 60F isodrosotherm if you want:

        #td60_c = ax1.contour(x,y,td2m,colors='darkgreen', alpha = 0.8, levels = [60])
        #td60_c.clabel(fontsize=12,colors='darkgreen',fmt='Td=60F',rightside_up=True,use_clabeltext=True)

        #Plot MSLP contours
        h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1080,2),linewidths=1,alpha=0.7)
        #Add label if you want:
        #h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

        #Define levels and colors for plotting precip types
        ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
        qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00','#ff9933','#ff6600','#ff5050','#ff0066']
        qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
        qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
        qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']

        capelevs = [100, 250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2500, 3000]

        #Plot surface-based CAPE
        capep = ax1.contourf(x, y, cape, levels=capelevs, alpha = 0.6, cmap='RdPu',extend='max')#['#0099ff00', '#4066ffb3', '#8066ff8c', '#BF66ff66','#8cff66','#b3ff66','#d9ff66','#ffff66','#ffd966','#ffcc66','#ffb366','#ff8c66','#ff6666','#ff668c','#ff66b3','#ff66d9','#ff66ff'])
        #Plot composite reflectivity
        ra = ax1.contourf(x,y,reflectivity,colors=qr_cols,levels=ref_levs,alpha=0.4,extend='max')

        #Plot 10m wind barbs
        ax1.barbs(x[wind_slice],y[wind_slice],u10m[wind_slice,wind_slice],v10m[wind_slice,wind_slice], length=6,color='gray')

        #Plot soundings
        smap.plot_soundings(fig,ax1,prs_temps,prs_relh,sfc_pressure,centerlat,centerlon,domainsize,model,cape=cape_bool,wetbulb=wetbulb_bool)

        #Add colorbars for CAPE and reflectivity
        spt.addcapecolorbar(ax1,fig,capep,capelevs)
        spt.addrefcolorbar(ax1,fig,ra,ref_levs)

        #Set plot extent and titles
        ax1.set_extent((west,east,south,north))
        ax1.set_title('Composite Reflectivity, SBCAPE, MSLP, 10m Wind, and Selected Soundings',fontsize=12)
        ax1.set_title('Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
        ax1.set_title(model+' Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')

        #Save plot
        plt.savefig(output_dir+'/'+model+'/'+domainname+'_severe_soundingmap_v22_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
        plt.close()
        plt.clf()
