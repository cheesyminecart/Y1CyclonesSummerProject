import tropycal as tp 
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 
import datetime as dt
import tropycal.tracks as tracks

# this tropycal module can be used for basin-wide tracks plotting, and could do nicely for a global warming investigation
# initialise data

wpac = tracks.TrackDataset(basin='west_pacific',source='ibtracs',ibtracs_mode='jtwc_neumann', include_btk=False)
wpac.to_dataframe()

# testing track resolution 

# wpac.plot_storm(('bolaven', 2023), domain='dynamic_tropical',prop={'dots':False,'linecolor':'category','linewidth':3.0})
# plt.show()

# ACE investigation, these are likely to be the bounds for our final dataset

# wpac.climatology(climo_bounds=(1981,2024)) 

# sample climatology function

wpac.ace_climo(plot_year=2023, compare_years=2017, climo_bounds=(1990, 2023))
plt.show()