import numpy as np
import pandas as pd
import tropycal.tracks as tracks
import datetime as dt

ibtracs = tracks.TrackDataset(
    basin='west_pacific',
    source='ibtracs',
    # ibtracs_url='/Users/bobo/Documents/GitHub/Y1CyclonesSummerProject/ibtracs.WP.list.v04r01.csv'
    ibtracs_url='https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r01/access/csv/ibtracs.WP.list.v04r01.csv'
)


storm = ibtracs.get_storm(('nepartak', 2021))
print(storm.to_dataframe())