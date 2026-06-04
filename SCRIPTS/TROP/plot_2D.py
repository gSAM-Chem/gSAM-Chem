import plotly.graph_objects as go
from netCDF4 import Dataset
import numpy as np

# Load the NetCDF file
dataset = Dataset('../../OUT_2D/TROP_4608x1152x72_SAM1MOM_4608_20200205_daily.2D_atm.nc')

# Extract the data
lons = dataset.variables['lon'][:]
lats = dataset.variables['lat'][:]
data = dataset.variables['LWNT'][:]

# Create a figure
fig = go.Figure(data=go.Heatmap(
                    z=data,
                    x=lons,
                    y=lats,
                    colorscale='Jet',
                ))

# Update layout to show the globe projection and coastlines
fig.update_geos(
    projection_type="natural earth",
    landcolor="white",
    oceancolor="MidnightBlue",
    showocean=True,
    lakecolor="LightBlue"
)

# Show the figure
fig.show()

