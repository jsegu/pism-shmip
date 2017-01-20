#!/usr/bin/env python2

import matplotlib.pyplot as plt
import netCDF4 as nc4

nc = nc4.Dataset('run-extra.nc')
x = nc.variables['x'][:]*1e-3
p = nc.variables['effbwp'][-1,:,:]
v = nc.variables['bwatvel[0]'][-1,:,:]
w = nc.variables['bwat'][-1,:,:]

# calc compute min, max and mean
pmin = p.min(axis=1)
pmax = p.max(axis=1)
pavg = p.mean(axis=1)
vavg = v.mean(axis=1)
wavg = w.mean(axis=1)
qavg = -wavg*vavg/(365.0*24*60*60)

# assert there are no longitudinal variations
assert (pmax-pmin).max() == 0.0

# init figure
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

# plot effective pressure
ax1.plot(x, pavg*1e-6)
ax1.set_ylim(0, 4.0)
ax1.set_ylabel('Effective pressure (MPa)')
ax1.grid()

# plot water velocity
ax2.plot(x, qavg*1e3)
ax2.set_xlim(0.0, 100.0)
ax2.set_ylim(0.0, 1.0)
ax2.set_xlabel('Distance from ice margin (km)')
ax2.set_ylabel('Water flux ($10^{-3}\,m^2\,s^{-1}$)')
ax2.grid()

# save
plt.savefig('plot.pdf')
