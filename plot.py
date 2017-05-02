#!/usr/bin/env python2

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

# use regular font for math text
plt.rc('mathtext', default='regular')
plt.rc('image', cmap='viridis')


def _envelope(x, z, ax=None, c=None, ls=None, axis=0):
    """Plot average and min to max envelope along given axis."""

    # get current axes if none
    ax = ax or plt.gca()

    # compute min, max and average profiles
    zmin = z.min(axis=axis)
    zmax = z.max(axis=axis)
    zavg = z.mean(axis=axis)

    # plot average and envelope
    l, = ax.plot(x, zavg, c=c, ls=ls)
    c = l.get_color()
    ax.fill_between(x, zmin, zmax, edgecolor='none', facecolor=c, alpha=0.25)


def plot_final(exp='a1'):
    """Plot y-min, y-max and y-avg final effective pressure and flux."""

    # experiment specific settins
    if exp[0] in ('e', 'f'):
        xmax = 6.0
    else:
        xmax = 100.0

    # open extra file
    print "Plotting experiment %s final stage..." % exp
    nc = nc4.Dataset('output/%s_extra.nc' % exp)

    # read config
    config = nc.variables['pism_config']
    g = config.standard_gravity
    c0 = config.till_cohesion
    phi = config.bootstrapping_tillphi_value_no_var
    rhoi = config.ice_density

    # read variables
    x = nc.variables['x'][:]*1e-3
    p = nc.variables['effbwp'][-1, :, :]*1e-6
    u = nc.variables['bwatvel[0]'][-1, :, :]
    v = nc.variables['bwatvel[1]'][-1, :, :]
    w = nc.variables['bwat'][-1, :, :]
    t = nc.variables['time'][-1]/(365.0*24*60*60)
    q = w*(u**2+v**2)**0.5/(365.0*24*60*60)*1e3

    # read till effective pressure
    if config.yield_stress_model == 'mohr_coulomb':
        n = nc.variables['tauc'][-1, :, :]*1e-6/np.tan(phi) - c0
    else:
        n = p * np.nan
    nc.close()

    # read overburden pressure from boot file
    if exp[0] == 'e':
        bfilename = 'input/boot_%s.nc' % exp[:2]
    elif exp[0] == 'f':
        bfilename = 'input/boot_e1.nc'
    else:
        bfilename = 'input/boot_sqrt.nc'
    nc = nc4.Dataset(bfilename)
    h = rhoi*g*nc.variables['thk'][0]*1e-6
    nc.close()

    # init figure
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # plot effective pressure
    ax1.set_title('Experiment %s, %.1f$\,a$' % (exp.upper(), t))
    _envelope(x, h, ax=ax1, c='0.5')
    _envelope(x, p, ax=ax1, axis=1)
    _envelope(x, n, ax=ax1, ls='--', axis=1)
    ax1.set_ylim(0.0, ax1.get_ylim()[1])
    ax1.set_ylabel('Effective pressure (MPa)')
    ax1.grid()

    # plot water velocity
    _envelope(x, q, ax=ax2, axis=1)
    ax2.set_xlim(0.0, xmax)
    ax2.set_ylim(0.0, ax2.get_ylim()[1])
    ax2.set_xlabel('Distance from ice margin (km)')
    ax2.set_ylabel('Water flux ($10^{-3}\,m^2\,s^{-1}$)')
    ax2.grid()

    # save
    plt.savefig('figures/final_%s' % exp)


def plot_timestamp(exp='a1'):
    """Plot real time versus model time."""

    # open extra file
    print "Plotting experiment %s time stamps..." % exp
    nc = nc4.Dataset('output/%s_extra.nc' % exp)
    t = nc.variables['time'][:]/(365.0*24*60*60)
    s = nc.variables['timestamp'][:]
    nc.close()

    # init figure
    fig, ax = plt.subplots(1, 1)

    # plot effective pressure
    ax.set_title('Experiment %s, %.2f$\,h$' % (exp.upper(), s[-1]))
    ax.plot(t, s)
    ax.set_xlabel('Model time (a)')
    ax.set_ylabel('Real time (h)')
    ax.grid()

    # save
    plt.savefig('figures/timestamp_%s' % exp)


def plot_transient(exp='a1'):
    """Plot time evolution of y-averaged effective pressure."""

    # experiment specific settins
    if exp[0] in ('e', 'f'):
        xmax = 6.0
    else:
        xmax = 100.0

    # open extra file
    print "Plotting experiment %s transient stage..." % exp
    nc = nc4.Dataset('output/%s_extra.nc' % exp)
    x = nc.variables['x'][:]*1e-3
    t = nc.variables['time'][:]/(365.0*24*60*60)
    p = nc.variables['effbwp'][:, :, :].mean(axis=2)*1e-6
    nc.close()

    # init figure
    fig, ax = plt.subplots(1, 1)

    # plot effective pressure
    ax.set_title('Experiment %s ' % exp.upper())
    im = ax.contourf(x, t, p)
    ax.set_xlabel('Distance from ice margin (km)')
    ax.set_ylabel('Time (a)')
    ax.set_xlim(0.0, xmax)
    ax.grid()

    # add colorbar
    cb = fig.colorbar(im)
    cb.set_label('Effective pressure (MPa)')

    # save
    plt.savefig('figures/transient_%s' % exp)


if __name__ == '__main__':
    """Main program, called at runtime."""

    # create inputs directory if absent
    if not os.path.exists('figures'):
        os.makedirs('figures')

    # parse arguments
    parser = argparse.ArgumentParser(description='Plot PISM SHMIP results.')
    parser.add_argument('-e', '--exp', type=str, default='a1',
                        help='experiment code (default: %(default)s)')
    args = parser.parse_args()

    # plot given experiment
    plot_final(args.exp)
    plot_timestamp(args.exp)
    plot_transient(args.exp)
