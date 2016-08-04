#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

from netCDF4 import Dataset, default_fillvals

fid = Dataset('../build_dbg/daily.nc', mode='r') # CLM one year test

# =========================
# Open and print dimensions
# =========================
ntime = len(fid.dimensions['time'])
print('ntime = '+str(ntime))
nlev = len(fid.dimensions['lev'])
print('nlev = '+str(nlev))

#lats_      = fid.variables['lat'][:]               # latitude array
#lons_      = fid.variables['lon'][:]               # longitude array
#times_ = fid.variables['time'][:] # time array

#time = fid.variables['time'][0:2*365]
NYEAR = 10

time = range(0,NYEAR*365)
zlev = fid.variables['zlev'][:]

#temp = fid.variables['T'][0:NYEAR*365,:]
#dens = fid.variables['rho'][0:NYEAR*365,:]
#print("Showing first "+str(NYEAR)+" years of simulation")

temp = fid.variables['T'][ntime-NYEAR*365:ntime,:]
dens = fid.variables['rho'][ntime-NYEAR*365:ntime,:]
print("Showing last "+str(NYEAR)+" years of simulation")

#cnt   = ax.contour(Tvec, Uvec, np.transpose(matSla), v, linewidths=0.5, colors='k')
fig, (ax1,ax2) = plt.subplots(1,2)
fig.suptitle('standalone firn model results')

X = time
Y = zlev

#my_cmap = plt.cm.get_cmap('inferno',21)
my_cmap = plt.cm.plasma


Z = np.transpose(temp)
im1         = ax1.contourf(X,Y,Z,11,cmap=my_cmap)
divider1    = make_axes_locatable(ax1)
cax1        = divider1.append_axes("bottom",size="10%", pad=0.25)
cbar1       = plt.colorbar(im1, cax=cax1, orientation="horizontal")

Z = np.transpose(dens)
im2         = ax2.contourf(X,Y,Z,11,cmap=my_cmap)
divider2    = make_axes_locatable(ax2)
cax2        = divider2.append_axes("bottom",size="10%", pad=0.25)
cbar2       = plt.colorbar(im2, cax=cax2, orientation="horizontal")

ax1.set_xlabel('day')
ax1.set_ylabel('depth')
ax1.invert_yaxis()
ax1.set_title('temperature')

ax2.set_xlabel('day')
ax2.set_ylabel('depth')
ax2.invert_yaxis()
ax2.set_title('density')

#plt.tight_layout()
plt.show()
