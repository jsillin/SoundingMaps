'''
This script hosts the function that plots soundings on maps.
'''

#### IMPORTS ####
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
import supplementary_tools as spt

def plot_soundings(fig,ax,temp,rh,sfc_pressure,centerlat,centerlon,domainsize,model,cape=False,wetbulb=False):
    """
    This function will plot a bunch of little soundings onto a matplotlib fig,ax.

    temp is an xarray dataarray with temperature data on pressure levels at least
    between 1000 and 300mb (you can change the ylimits for other datasets)

    rh is an xarray dataarray with temperature data on pressure levels at least
    between 1000 and 300mb (you can change )

    sfc_pressure is an xarray dataarray with surface pressure data (NOT MSLP!)

    centerlat and centerlon are the coordinates around which you want your map
    to be centered. both are floats or integers and are in degrees of latitude
    and degrees of longitude west (i.e. 70W would be input as positive 70 here)

    domainsize is a string either 'local' for ~WFO-size domains or 'regional' for
    NE/SE/Mid-Atlantic-size domains (12 deg lat by 15 deg lon). More will be added soon.

    model is a string that specifies which model is providing data for the plots.
    This determines a few things, most importantly longitude selections. Models
    currently supported are 'GFS','NAM',and 'RAP'

    cape is a boolean to indicate whether you want to overlay parcel paths and
    shade CAPE/CIN on soundings with >100 J/kg of CAPE (this value can be changed)

    wetbulb is a boolean to indicate whether you want to draw wet bulb profiles

    note that this function doesn't "return" anything but if you just call it and
    provide the right arguments, it works.

    for example:
        import soundingmaps as smap
        ...
        smap.plot_soundings(fig,ax1,data['temperature'],data['rh'],30.5,87.5,'local',cape=True)

    """
    r = 5
    if domainsize == "local":
        init_lat_delt = 1.625
        init_lon_delt = 0.45
        lat_delts = [0.2, 0.7, 1.2, 1.75, 2.25, 2.8]
        londelt = 0.76
        startlon = centerlon - 2 + 0.45

    elif domainsize == "regional":
        init_lat_delt = 6
        init_lon_delt = 1.6
        lat_delts = [0.6, 2.5, 4.5, 6.4, 8.4, 10.25]
        londelt = 2.9
        startlon = centerlon - 7.5 + 1.6

    # Lon adjustment for GFS because it's [0,360] not [-180,180]
    if model == 'GFS':
        startlon = 360-startlon

    # set lat/lon grid from which to pull data to plot soundings
    startlat = centerlat-init_lat_delt

    sound_lats = []
    sound_lons = []
    for i in range(0, 6):
        lats = startlat + lat_delts[i]
        sound_lats.append(lats)

    for i in range(0,r):
        if model == 'GFS':
            lons = startlon-(londelt*i)
        else:
            lons = -startlon-(londelt*i)
        sound_lons.append(lons)

    # this sets how high each row of soundings is on the plot
    plot_elevs=[0.2,0.3,0.4,0.5,0.6,0.7]

    # whole bunch of legend stuff
    dashed_red_line = lines.Line2D([], [], linestyle='solid', color='r', label='Temperature')
    dashed_purple_line = lines.Line2D([],[],linestyle='dashed',color='purple',label='0C Isotherm')
    dashed_green_line = lines.Line2D([], [], linestyle='solid', color='g', label='Dew Point')
    grey_line = lines.Line2D([], [], color='darkgray', label='MSLP (hPa)')
    blue_line = lines.Line2D([], [], color='b',label='Wet Bulb')
    pink_line = lines.Line2D([], [], color='fuchsia',label='Surface-Based Parcel Path')
    teal_line = lines.Line2D([], [], linestyle='dashed',color='teal',label='HGZ')
    green_dot = lines.Line2D([], [], marker='o', color='forestgreen',label='LCL')
    black_dot = lines.Line2D([], [], marker='o', color='k',label='Sounding Origin')

    red = mpatches.Patch(color='tab:red',label='CAPE')
    blue = mpatches.Patch(color='tab:blue',label='CIN')

    # do the plotting based on user inputs
    if cape and wetbulb is True:
            print('CAPE + Wetbulb')
            for i, plot_elev in enumerate(plot_elevs):
                soundlat = sound_lats[i]

                if k<2:
                    s=1
                else:
                    s=0

                for i in range(s,r):
                    levs_abv_ground = []
                    soundlon = sound_lons[i]
                    sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                    sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                    sound_pres = temp.lev
                    spres = sfc_pressure.interp(lat=soundlat,lon=soundlon)
                    sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
                    sound_wb = mpcalc.wet_bulb_temperature(sound_pres,sound_temps.data*units.degC,sound_dp)

                    #Only want data above the ground
                    abv_sfc_temp = spt.mask_below_terrain(spres,sound_temps,sound_pres)[0]
                    abv_sfc_dewp = spt.mask_below_terrain(spres,sound_dp,sound_pres)[0]
                    abv_sfc_wetb = spt.mask_below_terrain(spres,sound_wb,sound_pres)[0]
                    pres_abv_ground = spt.mask_below_terrain(spres,sound_temps,sound_pres)[1]

                    #sound_wb = sound_wb*units.degC
                    skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))

                    parcel_prof = mpcalc.parcel_profile(pres_abv_ground,abv_sfc_temp[0].data*units.degC,abv_sfc_dewp[0])
                    cape = mpcalc.cape_cin(pres_abv_ground,abv_sfc_temp.data*units.degC,abv_sfc_dewp,parcel_prof)
                    capeout = int(cape[0].m)
                    cinout = int(cape[1].m)

                    #skew.ax.axvspan(-30, -10, color='cyan', alpha=0.4)

                    skew.plot(pres_abv_ground,abv_sfc_wetb,'b',linewidth=2)
                    skew.plot(pres_abv_ground,abv_sfc_dewp,'g',linewidth=3)
                    skew.plot(pres_abv_ground,abv_sfc_temp,'r',linewidth=3)

                    if capeout >100:
                        # Shade areas of CAPE and CIN
                        print(pres_abv_ground)
                        print(abv_sfc_temp.data*units.degC)
                        print(parcel_prof)
                        skew.shade_cin(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof)
                        skew.shade_cape(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof)
                        skew.plot(pres_abv_ground,parcel_prof,color='fuchsia',linewidth=1)
                        lcl_pressure, lcl_temperature = mpcalc.lcl(pres_abv_ground[0], abv_sfc_temp.data[0]*units.degC, abv_sfc_dewp[0])
                        skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='forestgreen')
                        skew.ax.axvline(-30, color='teal', linestyle='--', linewidth=1)
                        skew.ax.axvline(-10, color='teal', linestyle='--', linewidth=1)
                    skew.plot(975,0, 'ko',markerfacecolor='k')

                    skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                    skew.ax.set_ylim((1000,300))
                    skew.ax.axis('off')

            leg = ax.legend(handles=[dashed_red_line,dashed_green_line,blue_line,dashed_purple_line,teal_line,green_dot,pink_line,red,blue,black_dot],title='Sounding Legend',loc=4,framealpha=1)

        else:
            print('CAPE no wetbulb')
            for k in range(len(plot_elevs)):
                soundlat = sound_lats[k]
                plot_elev = plot_elevs[k]

                if k==0:
                    s=1
                else:
                    s=0

                for i in range(s,r):
                    levs_abv_ground = []
                    soundlon = sound_lons[i]
                    sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                    sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                    sound_pres = temp.lev
                    spres = sfc_pressure.interp(lat=soundlat,lon=soundlon)
                    sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

                    abv_sfc_temp = spt.mask_below_terrain(spres,sound_temps,sound_pres)[0]
                    abv_sfc_dewp = spt.mask_below_terrain(spres,sound_dp,sound_pres)[0]
                    pres_abv_ground = spt.mask_below_terrain(spres,sound_temps,sound_pres)[1]

                    skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))

                    parcel_prof = mpcalc.parcel_profile(pres_abv_ground,abv_sfc_temp[0].data*units.degC,abv_sfc_dewp[0])
                    cape = mpcalc.cape_cin(pres_abv_ground,abv_sfc_temp.data*units.degC,abv_sfc_dewp,parcel_prof)
                    capeout = int(cape[0].m)
                    cinout = int(cape[1].m)

                    skew.plot(pres_abv_ground,abv_sfc_dewp,'g',linewidth=3)
                    skew.plot(pres_abv_ground,abv_sfc_temp,'r',linewidth=3)

                    if capeout >100:
                        # Shade areas of CAPE and CIN
                        skew.shade_cin(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof)
                        skew.shade_cape(pres_abv_ground, abv_sfc_temp.data*units.degC, parcel_prof)
                        skew.plot(pres_abv_ground,parcel_prof,color='fuchsia',linewidth=1)
                        print(abv_sfc_temp)
                        lcl_pressure, lcl_temperature = mpcalc.lcl(pres_abv_ground[0], abv_sfc_temp.data[0]*units.degC, abv_sfc_dewp[0])
                        skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='forestgreen')
                        skew.ax.axvline(-30, color='teal', linestyle='--', linewidth=1)
                        skew.ax.axvline(-10, color='teal', linestyle='--', linewidth=1)

                    skew.plot(975,0, 'ko',markerfacecolor='k')

                    skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                    skew.ax.set_ylim((1000,300))
                    skew.ax.axis('off')

            leg = ax.legend(handles=[dashed_red_line,dashed_green_line,dashed_purple_line,teal_line,green_dot,pink_line,red,blue,black_dot],title='Sounding Legend',loc=4,framealpha=1)

    else:
        if wetbulb==True:
            print('Wetbulb no CAPE')
            for k in range(len(plot_elevs)):
                soundlat = sound_lats[k]
                plot_elev = plot_elevs[k]

                if k==0:
                    s=1
                else:
                    s=0

                for i in range(s,r):
                    levs_abv_ground = []
                    soundlon = sound_lons[i]
                    sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                    sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                    sound_pres = temp.lev
                    spres = sfc_pressure.interp(lat=soundlat,lon=soundlon)

                    sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)

                    sound_wb = mpcalc.wet_bulb_temperature(sound_pres,sound_temps.data*units.degC,sound_dp)

                    abv_sfc_temp = spt.mask_below_terrain(spres,sound_temps,sound_pres)[0]
                    abv_sfc_dewp = spt.mask_below_terrain(spres,sound_dp,sound_pres)[0]
                    abv_sfc_wetb = spt.mask_below_terrain(spres,sound_wb,sound_pres)[0]
                    pres_abv_ground = spt.mask_below_terrain(spres,sound_temps,sound_pres)[1]

                    #sound_wb = sound_wb*units.degC
                    skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))

                    skew.plot(pres_abv_ground,abv_sfc_wetb,'b',linewidth=2)
                    skew.plot(pres_abv_ground,abv_sfc_dewp,'g',linewidth=3)
                    skew.plot(pres_abv_ground,abv_sfc_temp,'r',linewidth=3)

                    skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                    skew.ax.set_ylim((1000,300))
                    skew.ax.axis('off')
        else:
            print('No Wetbulb or CAPE')
            for k in range(len(plot_elevs)):
                soundlat = sound_lats[k]
                plot_elev = plot_elevs[k]

                if k==0:
                    s=1
                else:
                    s=0

                for i in range(s,r):
                    sound_pres = temp.lev
                    sound_temps = temp.interp(lat=soundlat,lon=soundlon)-273.15
                    sound_rh = rh.interp(lat=soundlat,lon=soundlon)
                    sound_dp = mpcalc.dewpoint_from_relative_humidity(sound_temps.data*units.degC,sound_rh.data*units.percent)
                    skew = SkewT(fig=fig,rect=(0.75-(0.15*i),plot_elev,.15,.1))
                    skew.plot(sound_pres,sound_dp,'g',linewidth=3)
                    skew.plot(sound_pres,sound_temps,'r',linewidth=3)
                    skew.plot(1000,0, 'ko',markerfacecolor='k')

                    skew.ax.axvline(0, color='purple', linestyle='--', linewidth=3)
                    skew.ax.set_ylim((1000,300))
                    skew.ax.axis('off')

            leg = ax.legend(handles=[dashed_red_line,dashed_green_line,blue_line,dashed_purple_line,black_dot],title='Sounding Legend',loc=4,framealpha=1)
