#!/usr/bin/env python2

import numpy as np
import netCDF4 as nc4


def init_pism_file(filename, x, y):
    """Init basic NetCDF file with x and y coords."""

    # open NetCDF file
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

    # close NetCDF file
    nc.close()


def make_boot_file(filename, x, y, b, h):
    """Make boot file with x and y coords and topg and thk variables."""

    # init NetCDF file
    print "Preparing boot file %s ..." % filename
    init_pism_file(filename, x, y)
    nc = nc4.Dataset(filename, 'a')

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

    # set ice surface temp
    var = nc.createVariable('ice_surface_temp', 'f4', ('y', 'x'))
    var[:] = 0.0*b + 260.0
    var.long_name = 'ice surface temperature for -surface given'
    var.units = 'Kelvin'

    # set surface mass balance
    var = nc.createVariable('climatic_mass_balance', 'f4', ('y', 'x'))
    var[:] = 0.0*b
    var.standard_name = 'land_ice_surface_specific_mass_balance_flux'
    var.long_name = 'climatic mass balance for -surface given'
    var.units = 'kg m-2 year-1'

    # set mask for prescribed sliding velocity
    var = nc.createVariable('bc_mask', 'f4', ('y', 'x'))
    var[:] = 0.0*b + 1.0
    var.long_name = 'mask prescribed sliding velocity'

    # set x-component of prescribed sliding velocity
    var = nc.createVariable('u_ssa_bc', 'f4', ('y', 'x'))
    var[:] = 0.0*b - 1e-6
    var.long_name = 'x-component of prescribed sliding velocity'
    var.units = 'm s-1'

    # set y-component of prescribed sliding velocity
    var = nc.createVariable('v_ssa_bc', 'f4', ('y', 'x'))
    var[:] = 0.0*b
    var.long_name = 'y-component of prescribed sliding velocity'
    var.units = 'm s-1'

    # close NetCDF file
    nc.close()


def make_boot_file_sqrt():
    """Make boot file for square root topography.

    To apply SHMIP-compliant boundary conditions we make two changes:

    * extend the domain by symetry in the x direction, and
    * add one ice-free grid cell immediately before x=0.
    """

    # prepare coordinates and topographies
    x = np.arange(-500.0, 200500.1, 500.0)
    y = np.arange(0.0, 20000.1, 500.0)
    xx, yy = np.meshgrid(x, y)
    xxsym = 100e3 - np.abs(xx-100e3)
    h = 1.0 + 6.0 * ((xxsym+5e3)**0.5-5e3**0.5)
    h[(xx<0.0)+(200000.0<xx)] = 0.0
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


def make_melt_file(filename, x, y, m):
    """Make basal melt input file with x and y coords and bmelt variable."""

    # init NetCDF file
    print "Preparing boot file %s ..." % filename
    init_pism_file(filename, x, y)
    nc = nc4.Dataset(filename, 'a')

    # set basal melt rate
    var = nc.createVariable('bmelt', 'f4', ('y', 'x'))
    var[:] = m
    var.standard_name = 'land_ice_basal_melt_rate'
    var.long_name = 'basal melt rate'
    var.units = 'm s-1'

    # close NetCDF file
    nc.close()


def make_melt_file_sqrt(filename, bgmelt=7.93e-11):
    """Make melt file for square root topography.

    To apply SHMIP-compliant boundary conditions we make two changes:

    * extend the domain by symetry in the x direction, and
    * add one ice-free grid cell immediately before x=0.
    """

    # prepare coordinates and topographies
    x = np.arange(-500.0, 200500.1, 500.0)
    y = np.arange(0.0, 20000.1, 500.0)
    xx, yy = np.meshgrid(x, y)
    xxsym = 100e3 - np.abs(xx-100e3)
    m = 0.0 * xx + bgmelt

    # make melt file
    make_melt_file(filename, x, y, m)


if __name__ == '__main__':
    """Main program, prepare all input files."""

    # prepare boot files
    make_boot_file_sqrt()
    make_boot_file_valley()

    # prepare melt files
    make_melt_file_sqrt('melt_a1.nc', bgmelt=7.93e-11)
    make_melt_file_sqrt('melt_a2.nc', bgmelt=1.59e-09)
    make_melt_file_sqrt('melt_a3.nc', bgmelt=5.79e-09)
    make_melt_file_sqrt('melt_a4.nc', bgmelt=2.5e-08)
    make_melt_file_sqrt('melt_a5.nc', bgmelt=4.5e-08)
    make_melt_file_sqrt('melt_a6.nc', bgmelt=5.79e-07)
