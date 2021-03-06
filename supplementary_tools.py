'''
This script contains tools for use in my various plotting routines, including
the soundingmaps.
'''

#### IMPORTS ####
from datetime import datetime
import datetime as dt
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.pyplot as plt

#### DATA FETCH HELPER FUNCTIONS ####
def get_init_time(model):
    '''
    This function will return date and run hour information to select the most
    current run of a given model.

    Input: model (string) currently supported 'HRRR','NAM','GFS','RTMA','RAP'

    Output: [mdate,init_hr] strings mdate=YYYYMMDD current run init_hr = HH.
    '''
    current_time = datetime.utcnow()
    year = current_time.year
    month = current_time.month
    day = current_time.day
    hour = current_time.hour

    if model=='HRRR':
        if hour <3:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<9:
            init_hour = '00'
        elif hour<14:
            init_hour = '06'
        elif hour<21:
            init_hour = '12'
        else:
            init_hour = '18'

    elif model=='NAM':
        if hour <4:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<10:
            init_hour = '00'
        elif hour<16:
            init_hour = '06'
        elif hour<22:
            init_hour = '12'
        else:
            init_hour = '18'

    elif model=='GFS':
        if hour <5:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            day = init_time.day
            month = init_time.month
            year = init_time.year
        elif hour<11:
            init_hour = '00'
        elif hour<17:
            init_hour = '06'
        elif hour<23:
            init_hour = '12'
        else:
            init_hour = '18'

    elif model=='RTMA':
        minute = current_time.minute
        if minute>50:
            init_hour = current_time.hour
            if float(init_hour) <10:
                init_hour = '0'+str(init_hour)
            else:
                init_hour = str(init_hour)
        else:
            time = current_time-dt.timedelta(hours=1)
            init_hour = time.hour
        if float(init_hour) <10:
            init_hour = '0'+str(init_hour)
        else:
            init_hour = str(init_hour)


    elif model=='RAP':
        minute = current_time.minute
        if minute<10:
            time = current_time-dt.timedelta(hours=2)
            init_hour = str(time.hour)
            if float(init_hour) <10:
                init_hour = '0'+str(init_hour)
            else:
                init_hour = str(init_hour)
        else:
            time = current_time-dt.timedelta(hours=1)
            init_hour = str(time.hour)
            if float(init_hour) <10:
                init_hour = '0'+str(init_hour)
            else:
                init_hour = str(init_hour)
    if month <10:
        month = '0'+str(month)
    else:
        month = str(month)

    if day <10:
        day = '0'+str(day)
    else:
        day = str(day)

    if hour <10:
        hour = '0'+str(hour)
    else:
        hour = str(hour)

    mdate = str(year)+month+day
    output = [mdate,init_hour]

    return output


def get_prev_init_time(model):
    '''
    This function will return date and run hour information for the previous
    forecast cycle of a given model. This is useful for analysis of model trends.

    Input: model (string) currently supported 'HRRR','NAM','GFS'

    Output: [mdate,init_hr] strings mdate=YYYYMMDD current run init_hr = HH.
    '''

    current_time = datetime.utcnow()
    year = current_time.year
    month = current_time.month
    day = current_time.day
    hour = current_time.hour

    if model=='HRRR':
        if hour <3:
            init_time = current_time-dt.timedelta(hours=3)
            init_hour = '18'
            prev_init_hour = '12'
            day = piday = init_time.day
            month = pimonth = init_time.month
            year = piyear= init_time.year
        elif hour<9:
            init_hour = '00'
            prev_init_hour = '18'
            prev_init_time = current_time-dt.timedelta(hours=9)
            piday = prev_init_time.day
            pimonth = prev_init_time.month
            piyear = prev_init_time.year
        elif hour<15:
            init_hour = '06'
            prev_init_hour = '00'
            piday = day
            pimonth = month
            piyear = year
        elif hour<21:
            init_hour = '12'
            prev_init_hour = '06'
            piday = day
            pimonth = month
            piyear = year
        else:
            init_hour = '18'
            prev_init_hour = '12'
            piday = day
            pimonth = month
            piyear = year

    elif model=='NAM':
        if hour <4:
            init_time = current_time-dt.timedelta(hours=4)
            init_hour = '18'
            prev_init_hour = '12'
            day = piday = init_time.day
            month = pimonth = init_time.month
            year = piyear= init_time.year
        elif hour<10:
            init_hour = '00'
            prev_init_hour = '18'
            prev_init_time = current_time-dt.timedelta(hours=10)
            piday = prev_init_time.day
            pimonth = prev_init_time.month
            piyear = prev_init_time.year
        elif hour<16:
            init_hour = '06'
            prev_init_hour = '00'
            piday = day
            pimonth = month
            piyear = year
        elif hour<22:
            init_hour = '12'
            prev_init_hour = '06'
            piday = day
            pimonth = month
            piyear = year
        else:
            init_hour = '18'
            prev_init_hour = '12'
            piday = day
            pimonth = month
            piyear = year
    elif model=='GFS':
        if hour <5:
            init_time = current_time-dt.timedelta(hours=5)
            init_hour = '18'
            prev_init_hour = '12'
            day = piday = init_time.day
            month = pimonth = init_time.month
            year = piyear= init_time.year
        elif hour<11:
            init_hour = '00'
            prev_init_hour = '18'
            prev_init_time = current_time-dt.timedelta(hours=11)
            piday = prev_init_time.day
            pimonth = prev_init_time.month
            piyear = prev_init_time.year
        elif hour<16:
            init_hour = '06'
            prev_init_hour = '00'
            piday = day
            pimonth = month
            piyear = year
        elif hour<22:
            init_hour = '12'
            prev_init_hour = '06'
            piday = day
            pimonth = month
            piyear = year
        else:
            init_hour = '18'
            prev_init_hour = '12'
            piday = day
            pimonth = month
            piyear = year

    if month <10:
        month = '0'+str(month)
    else:
        month = str(month)

    if day <10:
        day = '0'+str(day)
    else:
        day = str(day)

    if hour <10:
        hour = '0'+str(hour)
    else:
        hour = str(hour)

    if pimonth <10:
        pimonth = '0'+str(pimonth)
    else:
        pimonth = str(pimonth)

    if piday <10:
        piday = '0'+str(piday)
    else:
        piday = str(piday)

    mdate = str(piyear)+pimonth+piday
    output = [mdate,prev_init_hour]

    return output

def get_url(model):
    '''
    Return the NOMADS URL for a model of choice. Currently supported options are
    GFS, NAM, HRRR, RAP
    '''
    mdate = get_init_time(model)[0]
    init_hour = get_init_time(model)[1]
    if model == 'HRRR':
        url = 'http://nomads.ncep.noaa.gov:80/dods/hrrr/hrrr'+mdate+'/hrrr_sfc.t'+init_hour+'z'
    elif model == 'NAM':
        url = 'http://nomads.ncep.noaa.gov:80/dods/nam/nam'+mdate+'/nam_'+init_hour+'z'
    elif model == 'GFS':
        url = 'http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr/gfs'+mdate+'/gfs_0p25_1hr_'+init_hour+'z'
    elif model == 'RAP':
        url = 'http://nomads.ncep.noaa.gov:80/dods/rap/rap'+mdate+'/rap_'+init_hour+'z'

    return url

def get_num_timesteps(model):
    '''
    Return the number and width of time steps to query for a given model.
    Currently supported options are GFS, NAM, HRRR, RAP
    '''
    if model =='GFS':
        etime = 121
        delt = 1
    elif model == 'NAM':
        etime = 28
        delt = 3
    elif model == 'HRRR':
        etime = 49
        delt = 1
    elif model == 'RAP':
        etime = 37
        delt = 1

    return [etime,delt]

def get_varlist(model):
    '''
    Each model has slightly different variable names. This function will return
    a dictionary that renames the right variables to the right things depending
    on which model you want. Currently supported options are GFS, NAM, HRRR, RAP
    '''

    if model == 'RAP':
        vars = {
            'cfrzrsfc':'catice',
            'cicepsfc':'catsleet',
            'crainsfc':'catrain',
            'csnowsfc':'catsnow',
            'tmpprs': 'temperature',
            'mslmamsl':'mslp',
            'tmp2m':'sfc_temp',
            'dpt2m':'sfc_td',
            'refcclm':'radar',
            'rhprs':'rh',
            'capesfc':'cape',
            'ugrd10m':'u',
            'vgrd10m':'v',
            'pressfc':'spres'
        }

    elif model == 'HRRR':
        vars = {
            'cfrzrsfc':'catice',
            'cicepsfc':'catsleet',
            'crainsfc':'catrain',
            'csnowsfc':'catsnow',
            'tcdcclm':'tcc',
            'tmpprs': 'temperature',
            'ugrd10m': 'u',
            'vgrd10m': 'v',
            'mslmamsl':'mslp',
            'tmp2m':'sfc_temp',
            'dpt2m':'sfc_td',
            'refcclm':'radar',
            'apcpsfc':'qpf',
            'capesfc':'cape',
            'gustsfc':'sfcgust',
            'hcdchcll':'high_cloud',
            'mcdcmcll':'mid_cloud',
            'lcdclcll':'low_cloud',
            'vissfc':'sfcvis',
            'hgt263_k':'hgt_m10c',
            'hgt253_k':'hgt_m20c',
            'ltngclm':'lightning',
            'sbt124toa':'simsat',
            'hgt0c':'0chgt'
        }

    elif model == 'NAM':
        vars = {
            'cfrzrsfc':'catice',
            'cicepsfc':'catsleet',
            'crainsfc':'catrain',
            'csnowsfc':'catsnow',
            'tcdcclm':'tcc',
            'tmpprs': 'temperature',
            'ugrd10m': 'u',
            'vgrd10m': 'v',
            'hgtprs': 'height',
            'prmslmsl':'mslp',
            'tmp2m':'sfc_temp',
            'dpt2m':'sfc_td',
            'refcclm':'radar',
            'apcpsfc':'qpf',
            'rhprs':'rh',
            'capesfc':'cape',
            'pressfc':'spres'
        }

    elif model == 'GFS':
        vars = {
            'cfrzrsfc':'catice',
            'cicepsfc':'catsleet',
            'crainsfc':'catrain',
            'csnowsfc':'catsnow',
            'tcdcclm':'tcc',
            'tmpprs': 'temperature',
            'ugrd10m': 'u',
            'vgrd10m': 'v',
            'hgtprs': 'height',
            'prmslmsl':'mslp',
            'tmp2m':'sfc_temp',
            'dpt2m':'sfc_td',
            'refcclm':'radar',
            'apcpsfc':'qpf',
            'rhprs':'rh',
            'capesfc':'cape',
            'pressfc':'spres'

        }

    return vars

#### CALCULATION OF METEOROLOGICAL VARIABLES ####
def wet_bulb(temp,dewpoint):
    '''
    This uses the simple 1/3 rule to compute wet bulb temperatures from temp and
    dew point values/arrays. See Knox et. al (2017) in BAMS for more info about
    this approximation and when it is most reliable.

    Input: temp, dewpoint either values or arrays

    Output: wet_bulb either values or arrays depending on input
    '''
    tdd = temp-dewpoint
    wet_bulb = temp-((1/3)*tdd)
    return wet_bulb

def wetbulb_with_nan(pressure,temperature,dewpoint):
    '''
    This function uses the MetPy wet_bulb_temperature method to calculate the
    actual wet bulb temperature using pressure, temperature, and dew point info.

    Inputs: pressure, temperature, dewpoint pint arrays

    Output: wetbulb_full pint array

    This function was constructed using code graciously suggested by Jon Thielen
    '''
    nan_mask = np.isnan(pressure) | np.isnan(temperature) | np.isnan(dewpoint)
    idx = np.arange(pressure.size)[~nan_mask]
    wetbulb_valid_only = mpcalc.wet_bulb_temperature(pressure[idx], temperature[idx], dewpoint[idx])
    wetbulb_full = np.full(pressure.size, np.nan) * wetbulb_valid_only.units
    wetbulb_full[idx] = wetbulb_valid_only

    return wetbulb_full

def fram(ice,wet_bulb,velocity):
    '''
    This function computes ice accretion values using the Freezing Rain Accumulation
    Model method outlined in Sanders and Barjenbruch (2016) in WAF.

    Inputs: ice, wet_bulb, velocity which are arrays containing QPF falling as
    ZR, wet bulb temperature, and wind speed information. Units are inches per hour,
    degrees celsius, and knots respectively.

    Output: ice accretion array in units of inches.
    '''
    ilr_p = ice
    ilr_t = (-0.0071*(wet_bulb**3))-(0.039*(wet_bulb**2))-(0.3904*wet_bulb)+0.5545
    ilr_v = (0.0014*(velocity**2))+(0.0027*velocity)+0.7574

    cond_1 = np.ma.masked_where(wet_bulb>-0.35,ice)
    cond_2 = np.ma.masked_where((wet_bulb<-0.35) & (velocity>12.),ice)
    cond_3 = np.ma.masked_where((wet_bulb<-0.35) & (velocity<=12.),ice)

    cond_1 = cond_1.filled(0)
    cond_2 = cond_2.filled(0)
    cond_3 = cond_3.filled(0)

    ilr_1 = (0.7*ilr_p)+(0.29*ilr_t)+(0.01*ilr_v)
    ilr_2 = (0.73*ilr_p)+(0.01*ilr_t)+(0.26*ilr_v)
    ilr_3 = (0.79*ilr_p)+(0.2*ilr_t)+(0.01*ilr_v)

    accretion_1 = cond_1*ilr_1
    accretion_2 = cond_2*ilr_2
    accretion_3 = cond_3*ilr_3

    total_accretion=accretion_1+accretion_2+accretion_3
    return total_accretion

#### FIGURE TOOLS ####
def addcapecolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for CAPE

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0
    bottom = 0.17
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Surface-Based CAPE (J/kg)', size=8)  # MODIFY THIS for other fields!!

def addrefcolorbar(ax,fig,im,clevs):
    '''
    This function adds a new colorbar on its own axes for reflectivity

    Inputs: ax, fig are matplotlib axis/figure objects, im is the contourf object,
    and clevs is the contour levels used in the contourf plot

    Outputs: just call it and it'll put the colorbar in the right place

    This code was adapted from Dr. Kim Wood's Community Tools repo
    '''
    axes_bbox = ax.get_position()
    left = axes_bbox.x0 + 0.39
    bottom = 0.17
    width = 0.38
    height = 0.01
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=clevs, orientation='horizontal')
    #cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Composite Reflectivity (dBZ)', size=8)  # MODIFY THIS for other fields!!

#### MISC ####
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

def mask_below_terrain(spres,data,levs):
    '''
    Given a surface pressure, return data only below that pressure (above ground).

    Needs spres, a surface pressure (float)
    Needs data, a pint quantity array of temps/dew point/rh/whatever
    Needs levs, a pint quantity array of pressures
    '''
    above_ground = []
    for i in range(len(levs)):
        diff = levs[i]-spres
        if diff <0:
            above_ground.append(levs[i])
    pres_abv_ground = above_ground*units.hPa
    num_points_abv_ground = len(above_ground)
    data_abv_ground = data[-num_points_abv_ground:]
    return [data_abv_ground,pres_abv_ground]

def get_windslice(model,domainsize):
    '''
    Given a model and domainsize, return the right wind slice to make the barbs
    look good.

    Inputs: model, domainsize (strings). Currently supported models are 'GFS',
    'NAM', and 'RAP'. Currently supported domainsizes are 'regional' and 'local'

    Output: slice object wind_slice
    '''
    if model == 'GFS':
        if domainsize == 'regional':
            wind_slice = slice(2,-2,2)
        elif domainsize == 'local':
            wind_slice = slice(1,-1,1)
        else:
            print("Invalid domainsize String. Needs to be 'regional' or 'local'")
    elif model == 'NAM':
        if domainsize == 'regional':
            wind_slice = slice(6,-6,6)
        elif domainsize == 'local':
            wind_slice = slice(3,-3,3)
        else:
            print("Invalid domainsize String. Needs to be 'regional' or 'local'")
    elif model == 'RAP':
        if domainsize == 'regional':
            wind_slice = slice(12,-12,12)
        elif domainsize == 'local':
            wind_slice = slice(8,-8,8)
        else:
            print("Invalid domainsize String. Needs to be 'regional' or 'local'")
    else:
        print("Invalid model String. Needs to be 'GFS','NAM',or 'RAP'")

    return wind_slice
