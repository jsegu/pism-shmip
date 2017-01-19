#!/usr/bin/env python2

import matplotlib.pyplot as plt
import netCDF4 as nc4

nc = nc4.Dataset('run-extra.nc')
x = nc.variables['x'][:]*1e-3
p = nc.variables['effbwp'][-1,:,:]
v = nc.variables['bwatvel[0]'][-1,:,:]

# calc compute min, max and mean
pmin = p.min(axis=1)
pmax = p.max(axis=1)
pavg = p.mean(axis=1)
vavg = v.mean(axis=1)

# assert there are no longitudinal variations
assert (pmax-pmin).max() == 0.0

# plot effective pressure
plt.subplot(211)
plt.plot(x, pavg*1e-6)
plt.ylim(0, 4.0)
plt.ylabel('Water effective pressure (MPa)')

# plot water velocity
plt.subplot(212)
plt.plot(x, vavg*1e-6)
plt.ylim(-5.0, 5.0)
plt.ylabel('Water velocity ($10^6\,m\,a^{-1}$)')
plt.xlabel('Distance from ice margin (km)')

# save
plt.savefig('plot.pdf')
