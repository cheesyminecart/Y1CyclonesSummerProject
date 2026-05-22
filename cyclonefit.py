import numpy as np
import scipy as sc 
import xarray as xr
import matplotlib.pyplot as plt  # plotting library
# fitting file, with coordinate transformations, using testfile tammy_peakintensity.nc
# lat/long coordinate transformations per value, in kilometers 

latitude = 27.83 
longitude = 27.83 * np.cos(26.1) # using the relative latitude! 

# opening the file and reading it! 
tammy = xr.open_dataset(r'/path/to/file/Y1CyclonesSummerProject/data/tammy_peakintensity.nc')
tammy

