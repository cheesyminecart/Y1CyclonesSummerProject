#Comparing Hurricane Floyd to Hurricane Harvey using HURDAT2 

from tropycal import tracks, rain
import numpy as np
import matplotlib.pyplot as plt

# load HURDAT2 Dataset (North Atlantic)
basin = tracks.TrackDataset(basin='north_atlantic',include_btk=False)

rain_obj = rain.RainDataset()


storm_before = basin.get_storm(('floyd', 1999))   
storm_now    = basin.get_storm(('harvey', 2017))  

rain_before = rain_obj.get_storm_rainfall(storm_before)
rain_now    = rain_obj.get_storm_rainfall(storm_now)

# print maximum rainfall and other stats for each storm
def print_max_rain(df, name):
    row = df.loc[df['Total'] == np.nanmax(df['Total'])]
    print(f"--- {name} ---")
    print(f"Max Rainfall = {row['Total'].values[0]} Inches")
    print(f"Location     = {row['Station'].values[0]}")
    print(f"Latitude     = {row['Lat'].values[0]}")
    print(f"Longitude    = {row['Lon'].values[0]}")
    print()

print_max_rain(rain_before, "Storm BEFORE (Floyd 1999)")
print_max_rain(rain_now,    "Storm NOW (Harvey 2017)")

# plot gridded rainfall for both storms
grid_before = rain_obj.interpolate_to_grid(storm_before, return_xarray=True)
grid_now    = rain_obj.interpolate_to_grid(storm_now,    return_xarray=True)

levels = [1,2,4,8,12,16,20,30,40,50,60]

# Floyd (1999)
rain_obj.plot_rain_grid(
    storm_before,
    grid_before,
    levels,
    domain={'s':30, 'n':42, 'w':-85, 'e':-70}
)
plt.title("Rainfall Grid — Floyd (1999)")

# Harvey (2017)
rain_obj.plot_rain_grid(
    storm_now,
    grid_now,
    levels,
    domain={'s':26, 'n':39, 'w':-103, 'e':-83}
)
plt.title("Rainfall Grid — Harvey (2017)")

plt.show()
