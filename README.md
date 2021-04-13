# SoundingMaps
 Overlay small soundings above a map to analyze spatial and temporal trends in thermodynamic profiles

## 0.2 Update
 Version 0.2 is here as of 3/6/21!
 NEW features: 
 * get_soundingmaps.py is now the main interface through which you
 can make soundingmaps. Run in either interactive mode to set variables through
 the interpreter or tweak the variables in the script to pick your own region 
 and model of interest.
 * severe weather parameters:
   * CAPE/CIN shading
   * Parcel path plot 
   * HGZ outline for unstable soundings 
   * LCL marker for unstable soundings
 * sub-ground data now masked out 
 * option to plot wet bulb profiles
 * black dot marks origin of sounding data

This code now relies on MetPy 1.0, so make sure you're working with the most up-to-date version

 This is very much a work in progress! The plotting routines are a bit clumsy
 and it's not as easy as it should be to change domain locations/sizes.

 I'm posting this so others in the community can adapt it for their own ideas,
 help make it better, and also so I can learn more about how the development
 of open-source code actually works. Any and all tips/advice about how to work
 best in github, good coding practices, how to make this script or others like
 it more efficient, etc. are very much welcomed!

 Enjoy and happy coding.
 -Jack
