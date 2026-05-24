# simple algo to get storm coordinates, reducing data file size

import numpy as np
import pandas as pd
import tropycal.tracks as tracks
import datetime as dt

hurdat = tracks.TrackDataset(basin='north_atlantic',source='hurdat')

storm = hurdat.get_storm(('tammy',2023))
print(storm.to_dataframe())

