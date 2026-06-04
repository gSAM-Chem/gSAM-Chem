#!/usr/bin/env python

# --------------------------------------------------
# Non-interactive plotting (required on NCAR)
# --------------------------------------------------
import matplotlib
matplotlib.use("Agg")

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# --------------------------------------------------
# Input / output
# --------------------------------------------------
ncfile  = "out.nc"
outfile = "truck_W_QV_time_height.png"

# --------------------------------------------------
# Open NetCDF
# --------------------------------------------------
nc = Dataset(ncfile, "r")

W  = nc.variables["W"][:]         # (time,z) or (z,time)
QV = nc.variables["QV"][:]        # same shape
time = nc.variables["time"][:]    # irregular, 1D
z    = nc.variables["z"][:] * 0.001  # km

dims = nc.variables["W"].dimensions

# Ensure arrays are (z, time)
if dims[0] == "time":
    W  = W.T
    QV = QV.T

# Convert masked → NaN
W  = np.ma.filled(W,  np.nan)
QV = np.ma.filled(QV, np.nan)

# --------------------------------------------------
# Plot
# --------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    nrows=2, ncols=1,
    figsize=(10, 5),
    sharex=True
)

# ---- W panel ----
levels_W = np.arange(-5.0, 5.1, 0.1)

cf1 = ax1.contourf(
    time, z, W,
    levels=levels_W,
    cmap="RdBu_r",
    extend="both"
)

cbar1 = fig.colorbar(cf1, ax=ax1)
cbar1.set_label("W (m s$^{-1}$)")

ax1.set_ylabel("Height (km)")
ax1.set_title("Vertical velocity (W)")

# ---- QV panel ----
levels_QV = np.arange(0.0, 11.1, 0.1)

cf2 = ax2.contourf(
    time, z, QV,
    levels=levels_QV,
    cmap="viridis",
    extend="max"
)

cbar2 = fig.colorbar(cf2, ax=ax2)
cbar2.set_label("QV (g kg$^{-1}$)")

ax2.set_ylabel("Height (km)")
ax2.set_xlabel("Time (UTC)")
ax2.set_title("Specific humidity (QV)")

# ---- Time formatter: HH:MM ----
def hour_formatter(x, pos):
    h = int(x)
    m = int(round((x - h) * 60))
    return f"{h:02d}:{m:02d}"

ax2.xaxis.set_major_formatter(mticker.FuncFormatter(hour_formatter))

# --------------------------------------------------
# Save
# --------------------------------------------------
plt.tight_layout()
plt.savefig(outfile, dpi=200)
plt.close()

nc.close()

print(f"Wrote file: {outfile}")

