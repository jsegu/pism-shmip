#!/usr/bin/env python2
# coding: utf-8

import os
import numpy as np
import netCDF4 as nc4

# parameters
auth = 'jseg'
inst = 'ETH Zürich'
name = 'Julien Seguinot'
exps = (['a%d' % i for i in range(1, 5)] +
        ['b%d' % i for i in range(1, 6)] +
        ['e%d' % i for i in range(1, 5)])


def copy_attributes(ivar, ovar):
    """Copy all attributes from ivar to ovar."""
    for attname in ivar.ncattrs():
        setattr(ovar, attname, getattr(ivar, attname))


def postprocess(exp='a1'):
    """Convert PISM output to SHMIP conventions."""

    print "Postprocessing experiment %s..." % exp.upper()

    # boot filename
    if exp[0] == 'e':
        bfilename='input/boot_%s.nc' % exp
    else:
        bfilename='input/boot_sqrt.nc'

    # extra and output filenames
    efilename = 'output/%s_extra.nc' % exp
    ofilename = '%s/%s_%s.nc' % (auth, exp.upper(), auth)

    # open datasets
    bds = nc4.Dataset(bfilename, 'r')
    eds = nc4.Dataset(efilename, 'r')
    ods = nc4.Dataset(ofilename, 'w')

    # copy global attributes from extra file
    copy_attributes(eds, ods)

    # read coordinate variables from extra file
    x = eds.variables['x']
    y = eds.variables['y']
    t = eds.variables['time']
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # set additional global attributes
    ods.title = 'PISM experiment %s.' % exp.upper()
    ods.meshtype = 'structured'
    ods.dimension = '2D'
    ods.channels_on_edges = 'no'
    ods.institution = '%s, %s' % (name, inst)

    # create dimensions
    ods.createDimension('dim', 2)  # number of spatial dimensions
    ods.createDimension('xdim', len(x)*len(y))  # regular grid
    ods.createDimension('xdims', len(x)*len(y))  # staggered grid
    ods.createDimension('time', None)

    # create new coordinate variables
    ovar = ods.createVariable('xy', x.dtype, ('dim', 'xdim'))
    ovar[:] = np.meshgrid(x, y)
    ovar = ods.createVariable('xys', x.dtype, ('dim', 'xdims'))
    ovar[:] = np.meshgrid(x+dx/2, y+dy/2)
    ovar = ods.createVariable('Dxy', x.dtype, ('dim'))
    ovar[:] = (dx, dy)
    ovar = ods.createVariable('time', t.dtype, ('time'))
    ovar[:] = t[:]

    # copy boot bedrock topography
    bvar = bds.variables['topg']
    ovar = ods.createVariable('B', bvar.dtype, ('xdim'))
    ovar[:] = bvar[:]
    copy_attributes(bvar, ovar)

    # copy boot ice thickness
    bvar = bds.variables['thk']
    ovar = ods.createVariable('H', bvar.dtype, ('xdim'))
    ovar[:] = bvar[:]
    copy_attributes(bvar, ovar)

    # copy effective pressure
    evar = eds.variables['effbwp']
    ovar = ods.createVariable('N', evar.dtype, ('time', 'xdim'))
    ovar[:] = evar[:]
    copy_attributes(evar, ovar)

    # close datasets
    bds.close()
    eds.close()
    ods.close()


if __name__ == '__main__':
    """Main program, prepare all input files."""

    # create inputs directory if absent
    if not os.path.exists(auth):
        os.makedirs(auth)

    # postprocess experiments completed
    for exp in exps:
        postprocess(exp)
