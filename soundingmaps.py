import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import metpy
from datetime import datetime
import datetime as dt
from metpy.units import units
from metpy.plots import SkewT
import matplotlib.patches as mpatches
import matplotlib.lines as lines
import metpy.calc as mpcalc

def plot_soundings(fig,ax,temp,rh,centerlat,centerlon,domainsize,cape):
    '''
    This function will plot a bunch of little soundings onto a matplotlib fig,ax.

    temp is an xarray dataarray with temperature data on pressure levels at least
    between 1000 and 300mb (you can change the ylimits for other datasets)

    rh is an xarray dataarray with temperature data on pressure levels at least
    between 1000 and 300mb (you can change )

    centerlat and centerlon are the coordinates around which you want your map
    to be centered. both are floats or integers and are in degrees of latitude
    and degrees of longitude west (i.e. 70W would be input as positive 70 here)

    domainsize is a string either 'local' for ~WFO-size domains, 'regional' for
    NE/SE/Mid-Atlantic-size domains, or 'overview' to see about half the CONUS.

    cape is a boolean to indicate whether you want to overlay parcel paths and
    shade CAPE/CIN on soundings with >100 J/kg of CAPE (this value can be changed)

    note that this function doesn't "return" anything but if you just call it and
    provide the right arguments, it works.

    for example:
        import soundingmaps as smap
        ...
        smap.plot_soundings(fig,ax1,data['temperature'],data['rh'],30.5,87.5,'local',cape=True)

    '''
    r=5
    if domainsize=='local':
        init_lat_delt = 1.625
        init_lon_delt = 0.45
        lat_delts = [.2,.7,1.2,1.75,2.25,2.8]
        londelt = 0.76
    else:
        print('domainsize not supported yet')

    startlat = centerlat-init_lat_delt
    startlon = centerlon-2+0.45

    sound_lats=[]
    sound_lons=[]
    for i in range(0,6):
        lats = startlat+lat_delts[i]
        sound_lats.append(lats)

    for i in range(0,r):
        lons = -startlon-(londelt*i)
        sound_lons.append(lons)

    plot_elevs=[0.2,0.3,0.4,0.5,0.6,0.7]

    dashed_red_line = lines.Line2D([], [], linestyle='solid', color='r', label='Temperature')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='0C Isotherm')
    dashed_green_line = lines.Line2D([], [], linestyle='solid', color='g', label='Dew Point')
    grey_line = lines.Line2D([], [], color='darkgray', label='MSLP (hPa)')
    blue_line = lines.Line2D([], [], color='b',label='2m 0C Isotherm')
    pink_line = lines.Line2D([], [], color='fuchsia',label='Surface-Based Parcel Path')
    red = mpatches.Patch(color='tab:red',label='CAPE')
    blue = mpatches.Patch(color='tab:blue',label='CIN')

    if cape==True:
        for k in range(len(plot_elevs)):
            soundlat = sound_lats[k]
            plot_elev = plot_elevs[k]

            if k==0:
                s=1
            else:
                s=0

            for i in range(s,r):
                sound_pres = temp.lev
                soundlon = -(startlon+(londelt*i))
                sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
                skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))

                parcel_prof = mpcalc.parcel_profile(sound_pres,sound_temps[0].data*units.degC,sound_dp[0])
                cape = mpcalc.cape_cin(sound_pres,sound_temps.data*units.degC,sound_dp,parcel_prof)
                capeout = int(cape[0].m)
                cinout = int(cape[1].m)

                skew.plot(sound_pres,sound_dp,'g',linewidth=3)
                skew.plot(sound_pres,sound_temps,'r',linewidth=3)

                if capeout >100:
                    # Shade areas of CAPE and CIN
                    print(sound_temps)
                    print(parcel_prof)
                    skew.shade_cin(sound_pres, sound_temps.data*units.degC, parcel_prof)
                    skew.shade_cape(sound_pres, sound_temps.data*units.degC, parcel_prof)
                    skew.plot(sound_pres,parcel_prof,color='fuchsia',linewidth=1)

                skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                skew.ax.set_ylim((1000,300))
                skew.ax.axis('off')

        leg = ax.legend(handles=[dashed_red_line,dashed_green_line,dashed_purple_line,pink_line,red,blue],title='Sounding Legend',loc=4,framealpha=1)

    else:
        for k in range(len(plot_elevs)):
            soundlat = sound_lats[k]
            plot_elev = plot_elevs[k]

            if k==0:
                s=1
            else:
                s=0

            for i in range(s,r):
                soundlon = -(startlon+(londelt*i))
                sound_pres = temp.lev
                sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
                skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))
                skew.plot(sound_pres,sound_dp,'g',linewidth=3)
                skew.plot(sound_pres,sound_temps,'r',linewidth=3)
                skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                skew.ax.set_ylim((1000,300))
                skew.ax.axis('off')
        leg = ax.legend(handles=[dashed_red_line,dashed_green_line,dashed_purple_line],title='Sounding Legend',loc=4,framealpha=1)
