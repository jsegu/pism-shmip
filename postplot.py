#!/usr/bin/env python2

import os
import argparse
import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np

# use regular font for math text
plt.rc('mathtext', default='regular')
plt.rc('image', cmap='viridis')


def plot_final(exp='a1'):
    """Plot effective pressure from postprocessed file."""

    # open postprocessed file
    pfilename = 'processed/%s_jseg.nc' % exp.upper()
    nc = nc4.Dataset(pfilename, mode='r')
    x, y = nc.variables['coords1'][:, :]
    p = nc.variables['N'][:, :]

    # plot profile along centerline
    center = np.where(y == 10000)
    plt.scatter(x[center], p[-1, center])

    # show result on screen
    plt.show()


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
