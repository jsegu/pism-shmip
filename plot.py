#!/usr/bin/env python2

import os
import argparse
import matplotlib.pyplot as plt
import netCDF4 as nc4

# use regular font for math text
plt.rc('mathtext', default='regular')


def plot_minmax(exp='a1'):
    """Plot min, max and avg effective pressure and flux."""

    # open extra file
    print "Plotting experiment %s ..." % exp
    nc = nc4.Dataset('output/%s_extra.nc' % exp)
    x = nc.variables['x'][:]*1e-3
    p = nc.variables['effbwp'][-1,:,:]*1e-6
    v = nc.variables['bwatvel[0]'][-1,:,:]
    w = nc.variables['bwat'][-1,:,:]
    t = nc.variables['time'][-1]/(365.0*24*60*60)
    q = -v*w/(365.0*24*60*60)*1e3
    nc.close()

    # compute min, max and mean
    pmin = p.min(axis=1)
    pmax = p.max(axis=1)
    pavg = p.mean(axis=1)
    qmin = q.min(axis=1)
    qmax = q.max(axis=1)
    qavg = q.mean(axis=1)

    # init figure
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

    # plot effective pressure
    ax1.set_title('Experiment %s, %.1f$\,a$' % (exp.upper(), t))
    ax1.fill_between(x, pmin, pmax, edgecolor='none', alpha=0.25)
    ax1.plot(x, pavg)
    ax1.set_ylim(0.0, ax1.get_ylim()[1])
    ax1.set_ylabel('Effective pressure (MPa)')
    ax1.grid()

    # plot water velocity
    ax2.fill_between(x, qmin, qmax, edgecolor='none', alpha=0.25)
    ax2.plot(x, qavg)
    ax2.set_xlim(0.0, 100.0)
    ax2.set_ylim(0.0, ax2.get_ylim()[1])
    ax2.set_xlabel('Distance from ice margin (km)')
    ax2.set_ylabel('Water flux ($10^{-3}\,m^2\,s^{-1}$)')
    ax2.grid()

    # save
    plt.savefig('figures/%s' % exp)


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
    plot_minmax(args.exp)
