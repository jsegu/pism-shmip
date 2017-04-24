#!/usr/bin/env python2

import os
import numpy as np
import netCDF4 as nc4

# day and year in seconds
day = 24.0 * 60.0 * 60.0
year = 365.0 * day


def get_time_coord(diurnal=False, seasonal=False):
    """Built time coordinate as constant, diurnal or seasonal."""

    # build time array
    if diurnal and seasonal:
        raise ValueError('Can not combine diurnal and seasonal cycles.')
    elif diurnal:
        # prepare periodic input as PISM can't do it with daily periods
        t = np.linspace(5*year, 5*year+73*day, 73*288+1)  # 5-min interval
    elif seasonal:
        t = np.linspace(0.0, year, 365+1)  # 1-day interval
    else:
        t = np.array([0.0])

    # build time bounds array
    if diurnal or seasonal:
        dt = t[1] - t[0]
        tb = np.array([t - dt/2, t + dt/2]).T
    else:
        tb = None

    # return time and bounds
    return t, tb


def get_moulins_melt(t, x, y, s, moulins_file, moulins_relamp=0.0):
    """Compute distributed moulins input from pointwise csv file."""

    # init melt and load moulins file data
    melt = np.zeros((len(t), len(y), len(x)))
    moulins = np.loadtxt(moulins_file, delimiter=',', ndmin=2)

    # apply distributed melt rates on nearest grid cells
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    for (km, xm, ym, qm) in moulins:
        jm = np.abs(y - ym).argmin()
        im = np.abs(x - xm).argmin()
        imsym = np.abs(x + xm - 200e3).argmin()
        meltseries = qm*(1 - moulins_relamp*np.sin(2*np.pi*t/day))
        meltseries = np.maximum(0, meltseries)/dx/dy
        melt[:, jm, im] += meltseries
        melt[:, jm, imsym] += meltseries

    # return spatio-temporal melt array
    return melt


def get_seasonal_melt(t, s, temp_offset=0.0, lapse_rate=-0.0075,
                      ddf=0.01/86400):
    """Compute seasonal melt using a simple degree-day model."""

    # compute seasonal cycle
    temp = -16.0*np.cos(2*np.pi*t/year) - 5.0 + temp_offset

    # apply temperature lapse rate
    temp = [s*lapse_rate] + temp[:, None, None]

    # compute melt by simple degree-day model
    melt = np.maximum(0, temp*ddf)

    # return spatio-temporal melt array
    return melt


def get_topographies(mode, para):
    """Get spatial coordinates and surface topography.

    To apply PISM-compliant boundary conditions we make two changes:

    * extend the domain by symetry in the x direction, and
    * add one ice-free grid cell immediately before x=0.
    """

    # square root topography mode
    if mode == 'sqrt':

        # prepare spatial coordinates
        xmax = 200e3
        ymax = 20e3
        dx = dy = 500.0
        x = np.arange(-dx, xmax + dx + 0.1, dx)
        y = np.arange(0.0, ymax + 0.1, dy)
        xx, yy = np.meshgrid(x, y)
        xxsym = xmax/2 - np.abs(xx-xmax/2)

        # prepare topographies
        s = 1.0 + 6.0 * ((xxsym+5e3)**0.5-5e3**0.5)
        s[(xx < 0.0)+(xmax < xx)] = 0.0
        b = 0.0 * s
        h = s

    # valley topography mode
    elif mode == 'valley':

        # prepare spatial coordinates
        xmax = 6000.0
        ymax = 550.0
        dx = dy = 20.0
        x = np.arange(-dx, xmax + dx + 0.1, dx)
        y = np.arange(-ymax, ymax + 0.1, dy)
        xx, yy = np.meshgrid(x, y)

        # prepare surface topography
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
        h[h < 0.0] = 0.0

    # return coordinates and surface topography
    return x, y, b, h, s


def init_pism_file(filename, x, y, t, tb=None):
    """Init basic NetCDF file with x and y coords."""

    # open NetCDF file
    nc = nc4.Dataset(filename, 'w')

    # define the dimensions
    nc.createDimension('time', None)
    nc.createDimension('x', len(x))
    nc.createDimension('y', len(y))

    # set time coordinate
    tvar = nc.createVariable('time', 'f8', ('time',))
    tvar[:] = t
    tvar.axis = 'T'
    tvar.long_name = 'time'
    tvar.standard_name = 'time'
    tvar.calendar = '365_day'
    tvar.units = 's'

    # set time bounds
    if tb is not None:
        nc.createDimension('nv', 2)
        nc.variables['time'].bounds = 'time_bounds'
        tbvar = nc.createVariable('time_bounds', 'f8', ('time', 'nv'))
        tbvar[:] = tb

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


def make_boot_file(filename, mode='sqrt', para=0.05):
    """Make boot file with x and y coords and topg and thk variables."""

    # get coordinates and topographies
    x, y, b, h, s = get_topographies(mode=mode, para=para)

    # init NetCDF file
    print "Preparing boot file %s ..." % filename
    init_pism_file(filename, x, y, 0)
    nc = nc4.Dataset(filename, 'a')

    # set bedrock surface elevation
    bvar = nc.createVariable('topg', 'f4', ('time', 'y', 'x'))
    bvar[0] = b
    bvar.long_name = 'bedrock surface elevation'
    bvar.standard_name = 'bedrock_altitude'
    bvar.units = 'm'

    # set land ice thickness
    hvar = nc.createVariable('thk', 'f4', ('time', 'y', 'x'))
    hvar[0] = h
    hvar.long_name = 'land ice thickness'
    hvar.standard_name = 'land_ice_thickness'
    hvar.units = 'm'

    # set ice surface temp
    var = nc.createVariable('ice_surface_temp', 'f4', ('time', 'y', 'x'))
    var[0] = 0.0*b + 260.0
    var.long_name = 'ice surface temperature for -surface given'
    var.units = 'Kelvin'

    # set surface mass balance
    var = nc.createVariable('climatic_mass_balance', 'f4', ('time', 'y', 'x'))
    var[0] = 0.0*b
    var.standard_name = 'land_ice_surface_specific_mass_balance_flux'
    var.long_name = 'climatic mass balance for -surface given'
    var.units = 'kg m-2 year-1'

    # set mask for prescribed sliding velocity
    var = nc.createVariable('bc_mask', 'f4', ('time', 'y', 'x'))
    var[0] = 0.0*b + 1.0
    var.long_name = 'mask prescribed sliding velocity'

    # set x-component of prescribed sliding velocity
    var = nc.createVariable('u_ssa_bc', 'f4', ('time', 'y', 'x'))
    var[0] = 0.0*b - 1e-6
    var.long_name = 'x-component of prescribed sliding velocity'
    var.units = 'm s-1'

    # set y-component of prescribed sliding velocity
    var = nc.createVariable('v_ssa_bc', 'f4', ('time', 'y', 'x'))
    var[0] = 0.0*b
    var.long_name = 'y-component of prescribed sliding velocity'
    var.units = 'm s-1'

    # close NetCDF file
    nc.close()


def make_melt_file(filename, mode='sqrt', para=0.05, bgmelt=0.0,
                   moulins_file=None, moulins_relamp=0.0, temp_offset=None):
    """Make basal melt input file with x and y coords and bmelt variable."""

    # get time coordinate depending on options
    diurnal = (moulins_relamp != 0.0)
    seasonal = (temp_offset is not None)
    t, tb = get_time_coord(diurnal=diurnal, seasonal=seasonal)

    # get spatial coordinates and topographies
    x, y, b, h, s = get_topographies(mode=mode, para=para)

    # compute melt
    m = np.ones((len(t), len(y), len(x))) * bgmelt
    if seasonal:
        m += get_seasonal_melt(t, s, temp_offset=temp_offset)
    if moulins_file:
        m += get_moulins_melt(t, x, y, s, moulins_file, moulins_relamp)

    # init NetCDF file
    print "Preparing boot file %s ..." % filename
    init_pism_file(filename, x, y, t, tb)
    nc = nc4.Dataset(filename, 'a')

    # set time-dependent basal water input
    var = nc.createVariable('inputtobed', 'f4', ('time', 'y', 'x'))
    var[:] = m
    var.long_name = 'time-dependent basal water input'
    var.units = 'm s-1'

    # close NetCDF file
    nc.close()


if __name__ == '__main__':
    """Main program, prepare all input files."""

    # create inputs directory if absent
    if not os.path.exists('input'):
        os.makedirs('input')

    # prepare boot files
    make_boot_file('input/boot_sqrt.nc', mode='sqrt')
    make_boot_file('input/boot_e1.nc', mode='valley', para=0.05)
    make_boot_file('input/boot_e2.nc', mode='valley', para=0.0)
    make_boot_file('input/boot_e3.nc', mode='valley', para=-0.1)
    make_boot_file('input/boot_e4.nc', mode='valley', para=-0.5)
    make_boot_file('input/boot_e5.nc', mode='valley', para=-0.7)

    # prepare melt files for exp. B1 to B5
    kwa = dict(mode='sqrt')
    make_melt_file('input/melt_b1.nc', moulins_file='moulins_b1.csv', **kwa)
    make_melt_file('input/melt_b2.nc', moulins_file='moulins_b2.csv', **kwa)
    make_melt_file('input/melt_b3.nc', moulins_file='moulins_b3.csv', **kwa)
    make_melt_file('input/melt_b4.nc', moulins_file='moulins_b4.csv', **kwa)
    make_melt_file('input/melt_b5.nc', moulins_file='moulins_b5.csv', **kwa)

    # prepare melt files for exp. C1 to C4
    kwa = dict(mode='sqrt', moulins_file='moulins_b5.csv')
    make_melt_file('input/melt_c1.nc', moulins_relamp=0.25, **kwa)
    make_melt_file('input/melt_c2.nc', moulins_relamp=0.5, **kwa)
    make_melt_file('input/melt_c3.nc', moulins_relamp=1.0, **kwa)
    make_melt_file('input/melt_c4.nc', moulins_relamp=2.0, **kwa)

    # prepare melt files for exp. D1 to D5
    make_melt_file('input/melt_d1.nc', mode='sqrt', temp_offset=-4.0)
    make_melt_file('input/melt_d2.nc', mode='sqrt', temp_offset=-2.0)
    make_melt_file('input/melt_d3.nc', mode='sqrt', temp_offset=0.0)
    make_melt_file('input/melt_d4.nc', mode='sqrt', temp_offset=2.0)
    make_melt_file('input/melt_d5.nc', mode='sqrt', temp_offset=4.0)

    # prepare melt files for exp. F1 to F5
    make_melt_file('input/melt_f1.nc', mode='valley', temp_offset=-6.0)
    make_melt_file('input/melt_f2.nc', mode='valley', temp_offset=-3.0)
    make_melt_file('input/melt_f3.nc', mode='valley', temp_offset=0.0)
    make_melt_file('input/melt_f4.nc', mode='valley', temp_offset=3.0)
    make_melt_file('input/melt_f5.nc', mode='valley', temp_offset=6.0)
