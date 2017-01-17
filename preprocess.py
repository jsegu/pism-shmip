#!/usr/bin/env python2

import numpy as np
import netCDF4 as nc4


def make_boot_file(filename, x, y, b, h):
    """Make boot file with x and y coords and topg and thk variables."""

    # open NetCDF file
    print "Preparing boot file %s ..." % filename
    nc = nc4.Dataset(filename, 'w')

    # define the dimensions
    nc.createDimension('x', len(x))
    nc.createDimension('y', len(y))

    # set projection x coordinate
    xvar = nc.createVariable('x', 'f8', ('x',))
    xvar[:] = x
    xvar.axis = 'X',
    xvar.long_name = 'x-coordinate in Cartesian system'
    xvar.standard_name = 'projection_x_coordinate'
    xvar.units = 'm'

    # set projection y coordinate
    yvar = nc.createVariable('y', 'f8', ('y',))
    yvar[:] = y
    yvar.axis = 'Y'
    yvar.long_name = 'y-coordinate in Cartesian system'
    yvar.standard_name = 'projection_y_coordinate'
    yvar.units = 'm'

    # set bedrock surface elevation
    bvar = nc.createVariable('topg', 'f4', ('y', 'x'))
    bvar[:] = b
    bvar.long_name = 'bedrock surface elevation'
    bvar.standard_name = 'bedrock_altitude'
    bvar.units = 'm'

    # set land ice thickness
    hvar = nc.createVariable('thk', 'f4', ('y', 'x'))
    hvar[:] = h
    hvar.long_name = 'land ice thickness'
    hvar.standard_name = 'land_ice_thickness'
    hvar.units = 'm'

    # close NetCDF file
    nc.close()


def make_boot_file_sqrt():
    """Make boot file for square root topography."""

    # prepare coordinates and topographies
    x = np.arange(0.0, 100000.1, 500.0)
    y = np.arange(0.0, 20000.1, 500.0)
    xx, yy = np.meshgrid(x, y)
    h = 1.0 + 6.0 * ((xx+5e3)**0.5-5e3**0.5)
    b = 0.0 * h

    # make boot file
    make_boot_file('boot_sqrt.nc', x, y, b, h)


def make_boot_file_valley(para=0.05):
    """Make boot file for valley topography."""

    # prepare coordinates
    xmax = 6.0e3
    ymax = 550.0
    x = np.arange(0.0, xmax + 0.1, 20.0)
    y = np.arange(-ymax, ymax, 20.0)
    xx, yy = np.meshgrid(x, y)

    # surface topography
    s = 1.0 + 100.0*(xx/xmax+(xx+200.0)**0.25-200.0**0.25)
    smax = s[0, -1]

    # helper functions
    f_func = para*xx + (smax-para*xmax) / xmax**2 * xx**2
    f_benc = 0.05*xx + (smax-0.05*xmax) / xmax**2 * xx**2
    g_func = 0.5e-6 * abs(yy)**3
    h_func = (5 - 4.5*xx/xmax) * (s-f_func) / (s-f_benc+1e-12)

    # basal topography and thickness
    b = f_func + g_func * h_func
    h = s - b
    h[h<0.0] = 0.0

    # make boot file
    make_boot_file('boot_valley.nc', x, y, b, h)


if __name__ == '__main__':
    """Main program, prepare all input files."""

    # prepare boot files
    make_boot_file_sqrt()
    make_boot_file_valley()
