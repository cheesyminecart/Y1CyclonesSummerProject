# simple algo to get storm coordinates, reducing data file size

import numpy as np
import pandas as pd
import tropycal.tracks as tracks
import datetime as dt

<<<<<<< Updated upstream
hurdat = tracks.TrackDataset(basin='north_atlantic',source='hurdat')

storm = hurdat.get_storm(('isabel',2003))
=======
# ibtracs = tracks.TrackDataset(basin='west_pacific', source='ibtracs', ibtracs_url='path/to/ibtracs.WP.list.v04r01.csv')
ibtracs = tracks.TrackDataset(basin='western_pacific',source='ibtracs')

storm = ibtracs.get_storm(('nepartak',2021))
>>>>>>> Stashed changes
print(storm.to_dataframe())

