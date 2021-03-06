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
import matplotlib.patches as mpatches
import matplotlib.lines as lines

'''
This program produces soundings derived from GFS model data obtained via
the NOMADS openDAP functionality and overlays these soundings above a
map of precipitation type and MSLP to assist in the assessment of spatial
and temporal changes in thermodynamic and moisture profiles.

This code was originally written by Jack Sillin.
'''

# make unique directory to store output
def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    from errno import EEXIST
    from os import makedirs,path

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

#grabbing data from NOMADS
startTime=datetime.now()

year = startTime.year

if startTime.month <10:
    month = '0'+str(startTime.month)
else:
    month = str(startTime.month)

if startTime.day <10:
    day = '0'+str(startTime.day)
else:
    day = str(startTime.day)

if startTime.hour <10:
    hour = '0'+str(startTime.hour)
else:
    hour = str(startTime.hour)

mdate = str(year)+str(month)+str(day)

def get_init_hr(hour):
    if int(hour) <6:
        init_hour = '00'
    elif int(hour) <12:
        init_hour = '06'
    elif int(hour) <17:
        init_hour = '12'
    elif int(hour) <22:
        init_hour = '18'
    else:
        init_hour = '00'
    return(init_hour)

url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs'+mdate+'/gfs_0p25_1hr_'+get_init_hr(hour)+'z'
init_hour = get_init_hr(hour)

# Create new directory to store output
output_dir = str(year)+str(month)+str(day)+'_'+str(init_hour)+'00' #this string names the output directory
mkdir_p(output_dir)
mkdir_p(output_dir+'/GFS') #create subdirectory to store GFS output like this

#This actually opens the dataset from NOMADS and parses it with MetPy
ds = xr.open_dataset(url)
init_hr = dt.datetime(year,int(month),int(day),int(init_hour))
times = ds['tmp2m'].metpy.time #Pull out the time dimension
init_time = ds['time'][0]

#Subset the data to only work with certain lats and lons of interest
lats = np.arange(20,55,0.25)
lons = np.arange(240,310,0.25)

ds = ds.sel(lat = lats, lon = lons)

#Now loop through the 120 forecast hours to make the plots
for i in range(0,120):
    #Get the data for the forecast hour of interest
    data = ds.metpy.parse_cf()
    data = data.isel(time=i)

    #Rename variables to useful things
    data = data.rename({
        'cfrzrsfc':'catice',
        'cicepsfc':'catsleet',
        'crainsfc':'catrain',
        'csnowsfc':'catsnow',
        'tmpprs': 'temperature',
        'prmslmsl':'mslp',
        'tmp2m':'sfc_temp',
        'dpt2m':'sfc_td',
        'refcclm':'radar',
        'rhprs':'rh'
    })

    #Pull out the categorical precip type arrays
    catrain = data['catrain'].squeeze()
    catsnow = data['catsnow'].squeeze()
    catsleet = data['catsleet'].squeeze()
    catice = data['catice'].squeeze()

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
    vertical, = data['temperature'].metpy.coordinates('vertical')
    time = data['temperature'].metpy.time
    x, y = data['temperature'].metpy.coordinates('x', 'y')
    lat, lon = xr.broadcast(y, x)
    zH5_crs = data['temperature'].metpy.cartopy_crs

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

    #Process MSLP data
    mslp = data['mslp']/100.
    mslpc = mslp.squeeze()
    mslpc=ndimage.gaussian_filter(mslpc,sigma=1,order=0)

    #This creates a nice-looking datetime label
    dtfs = str(time.dt.strftime('%Y-%m-%d_%H%MZ').item())

    ########## SET UP FIGURE ##################################################
    fig = plt.figure(figsize=(15,15))
    ax1 = fig.add_subplot(111, projection = zH5_crs)

    ax1.coastlines(resolution='10m')
    ax1.add_feature(cfeature.BORDERS.with_scale('10m'))
    ax1.add_feature(cfeature.STATES.with_scale('10m'))

    ########## PLOTTING #######################################################
    #Plot 2m 32F isotherm
    tmp_2m32 = ax1.contour(x,y,t2m,colors='b', alpha = 0.8, levels = [32])

    #Plot labeled MSLP contours
    h_contour = ax1.contour(x, y, mslpc, colors='dimgray', levels=range(940,1080,4),linewidths=1,alpha=0.7)
    h_contour.clabel(fontsize=14, colors='dimgray', inline=1, inline_spacing=4, fmt='%i mb', rightside_up=True, use_clabeltext=True)

    #Define levels and colors for plotting precip types
    ref_levs = [1,5,10,15,20, 25, 30, 35, 40, 45, 50, 55, 60, 65]
    qr_cols = ['#cfffbf','#a7ff8a','#85ff5c','#60ff2b','#40ff00','#ffff00','#e6e600','#cccc00','#e4cc00']
    qs_cols = ['#b8ffff','#82ffff','#00ffff','#00cccc','#00adad','#007575','#0087f5','#0039f5','#1d00f5']
    qi_cols = ['#eeccff','#dd99ff','#cc66ff','#bb33ff','#aa00ff','#8800cc','#660099','#440066','#6600cc']
    qz_cols = ['#ff0066','#ff0080','#ff33cc','#ff00bf','#cc0099','#990073','#66004d','#b30000','#ff3333']

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

    #Set a title and extent for the map
    ax1.set_title('Precipitation Type and Selected Soundings',fontsize=16)
    ax1.set_title('\n Valid: '+time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='right')
    ax1.set_title('\n GFS Init: '+init_time.dt.strftime('%Y-%m-%d %H:%MZ').item(),fontsize=11,loc='left')
    ax1.set_extent((287.25, 291.25, 42.5, 45.75))#, crs = zH5_crs)    # Set a title and show the plot

    #################### SOUNDINGS ################################
    '''
    Ok this is the fun part. I admit this is probably the number one ugliest
    way of possibly doing this. So sorry for the frustration!

    Basically, we're going to run six for loops, one for each row of soundings.
    For each row, we're going to start a little west of the eastern plot bound
    and move west by some fixed increment that will depend on the size of your
    domain (the regional and local files will have different londelts).

    At each point of interest, we're going to interpolate the temp/rh fields
    at each pressure level and stuff that info into an array for plotting.
    Because we're using metpy, we also need to tag it with units. I'm told there
    are ways of doing this unit tagging that are 1000% better than what's shown
    here, but I always seem to find a way to break those. This way works, though
    it is very far from elegant.

    I've added a 0C isotherm with a purple dashed line since precip types are
    a pressing forecast concern at the time I wrote this (Feb 2021). I've also
    truncated the soundings at 300mb since stuff above that doesn't matter for
    precip types but you can change this with the ptop variable for severe
    applications.

    One other note of some importance is that because MetPy doesn't currently
    support data coordinates in the rect function, great care is required to
    make sure that the sounding data is actually coming from the bottom of your
    skew t plots (I set mine to come from the bottom of the purple line). I've
    put in a feature request for this, so hopefully it will be easier in the
    future.
    '''

    '''
    Local-specific stuff...
    So this sounding code is similar but a little different than the original.
    I've tried to make it a little easier to create different domains because
    manually lining up the locations of the plots with the locations from which
    data is being pulled is really time-consuming. To make your own regional
    domain without checking this, follow the following formulae:

    Step 1: decide on bounds for your area of interest. You have 4 degrees of
    longitude and 3.25 degrees of latitude to work with. Changing these numbers
    will require re-calibrating the plot locations! In set_extent above, put
    your bounds in.

    Step 2: change startlon below to be your eastern bound + 0.45 degrees. Leave
    londelt alone unless you want to recalibrate everything!

    Step 3: change startlat to be your southern bound (no addition needed)
    '''

    sound_pres = data.lev
    ptop=300
    startlon=69.2
    startlat=42.5
    londelt=0.76
    sound_lons = []
    sound_lats = []
    lat_delts = [.2,.7,1.2,1.75,2.25,2.8]
    r=5
    for i in range(0,r):
        lon = -startlon-(londelt*i)
        sound_lons.append(lon)

    for i in range(0,6):
        lat = startlat+lat_delts[i]
        sound_lats.append(lat)
    print(sound_lats)
    for i in range(1,r):
        soundlat = sound_lats[0]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.2,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    for i in range(0,r):
        soundlat = sound_lats[1]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.3,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    for i in range(0,r):
        soundlat = sound_lats[2]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.4,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    for i in range(0,r):
        soundlat = sound_lats[3]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.5,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    for i in range(0,r):
        soundlat = sound_lats[4]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.6,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    for i in range(0,r):
        soundlat = sound_lats[5]
        soundlon = 360-(startlon+(londelt*i))
        sound_temps = data['temperature'].interp(lat=soundlat,lon=soundlon)-273.15

        sound_rh = data['rh'].interp(lat=soundlat,lon=soundlon)
        sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

        skew = SkewT(fig=fig,rect=(0.75-(0.15*i),0.7,.15,.1))
        skew.plot(sound_pres,sound_dp,'g',linewidth=3)
        skew.plot(sound_pres,sound_temps,'r',linewidth=3)
        skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
        skew.ax.set_ylim((1000,ptop))
        skew.ax.axis('off')

    #uncomment the two lines below (rows,cols) to plot gridlines and check that
    #your sounding bases line up with the right lats/lons used to pull the data

    #rows = ax1.gridlines(ylocs=sound_lats,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)
    #cols = ax1.gridlines(xlocs=sound_lons,linewidth=2, linestyle='--', edgecolor='dimgrey',draw_labels=True)

    ########## LEGEND ############
    dashed_red_line = lines.Line2D([], [], linestyle='solid', color='r', label='Temperature')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='0C Isotherm')
    dashed_green_line = lines.Line2D([], [], linestyle='solid', color='g', label='Dew Point')
    grey_line = lines.Line2D([], [], color='darkgray', label='MSLP (hPa)')
    blue_line = lines.Line2D([], [], color='b',label='2m 0C Isotherm')
    leg = ax1.legend(handles=[dashed_red_line,dashed_green_line,dashed_purple_line],title='Sounding Legend',loc=4,framealpha=1)
    leg.set_zorder(100)

    ######## Save the plot
    plt.savefig(output_dir+'/GFS/gfs_hrly_pytpe_gyx_sound_'+dtfs+'.png',bbox_inches='tight',pad_inches=0.1)
    fcst_hr = str(0)
    plt.close()
    plt.clf()
