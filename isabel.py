
# done by JM
# for future reference please copy your code into separate blocks! 
import datetime as dt
from tropycal import tracks

# Load HURDAT2
basin = tracks.TrackDataset(basin='north_atlantic', include_btk=False)

# Get the Hurricane Isabel data
isabel = basin.get_storm(('isabel', 2003))


t = dt.datetime(2003, 9, 18, 18)

# Find index
idx = isabel.dict['time'].index(t)

# Find the coordinates of the centre of the cyclone
lat = isabel.dict['lat'][idx]
lon = isabel.dict['lon'][idx]
typ = isabel.dict['type'][idx]

print("Isabel centre at 18 Sep 2003 18:00 UTC:")
print("Latitude:", lat)
print("Longitude:", lon) 
print(typ)