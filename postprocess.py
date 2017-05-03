#!/usr/bin/env python2
# coding: utf-8

import os
import shutil
import numpy as np
import netCDF4 as nc4
import tarfile


# parameters
auth = 'jseg'
inst = 'ETH Zürich'
name = 'Julien Seguinot'
exps = (['a%d' % i for i in range(1, 7)] +
        ['b%d' % i for i in range(1, 6)] +
        ['c%d' % i for i in range(1, 5)] +
        ['d%d' % i for i in range(1, 6)] +
        ['e%d' % i for i in range(1, 6)] +
        ['f%d' % i for i in range(1, 6)])


# day and year in seconds
day = 24.0 * 60.0 * 60.0
year = 365.0 * day


def copy_attributes(ivar, ovar):
    """Copy all attributes from ivar to ovar. Rename long_name to pism_name."""
    for attname in ivar.ncattrs():
        newname = 'pism_name' if attname == 'long_name' else attname
        setattr(ovar, newname, getattr(ivar, attname))


def postprocess(exp='a1'):
    """Convert PISM output to SHMIP conventions."""

    # boot filename
    if exp[0] == 'e':
        bfilename = 'input/boot_%s.nc' % exp[:2]
    elif exp[0] == 'f':
        bfilename = 'input/boot_e1.nc'
    else:
        bfilename = 'input/boot_sqrt.nc'

    # extra and output filenames
    efilename = 'output/%s_extra.nc' % exp
    ofilename = 'processed/%s_%s.nc' % (exp[:2].upper() + exp[2:], auth)

    # check for file presence
    if not os.path.isfile(bfilename):
        print "Warning: could not find input for exp %s." % exp.upper()
        return
    elif not os.path.isfile(efilename):
        print "Warning: could not find output for exp %s." % exp.upper()
        return
    else:
        print "Postprocessing experiment %s..." % exp.upper()

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

    # prepare slicing on x coordinate
    if exp[0] in ('a', 'b', 'c', 'd'):
        xcond = (0e3 <= x[:]) * (x[:] <= 100e3)
    else:
        xcond = (0e3 <= x[:]) * (x[:] <= 6e3)

    # prepare slicing on time coordinate
    if exp[0] in ('c'):
        tcond = (t[:] >= t[-1] - day)
    elif exp[0] in ('d', 'f'):
        tcond = (t[:] >= t[-1] - year)
    else:
        tcond = (t[:] >= 0.0)

    # prepare sliced coordinates
    ys = y[:]
    xs = x[xcond]
    ts = t[tcond]

    # set additional global attributes
    ods.title = 'PISM experiment %s.' % exp.upper()
    ods.meshtype = 'structured'
    ods.dimension = '2D'
    ods.channels_on_edges = 'no'
    ods.institution = '%s, %s' % (name, inst)
    ods.references = ('http://pism-docs.org, '
                      'https://shmip.bitbucket.io, '
                      'https://github.com/jsegu/pism-shmip')

    # create dimensions
    ods.createDimension('dim', 2)  # number of spatial dimensions
    ods.createDimension('index1', len(xs)*len(ys))  # regular grid
    ods.createDimension('index2', len(xs)*len(ys))  # staggered grid
    ods.createDimension('time', None)

    # create SHMIP node (PISM cell center) coordinate variables
    ovar = ods.createVariable('coords1', x.dtype, ('dim', 'index1'))
    ovar[:] = np.meshgrid(ys, xs)[::-1]
    ovar.long_name = 'node coordinates'
    ovar.pism_name = 'cell center coordinate'
    ovar.units = x.units

    # create SHMIP cell (PISM staggered) coordinate variables
    # (I deduced the sign of x and y shifts by looking at model output)
    ovar = ods.createVariable('coords2', x.dtype, ('dim', 'index2'))
    ovar[:] = np.meshgrid(ys+dy/2, xs+dx/2)[::-1]
    ovar.long_name = 'cell midpoint coordinates'
    ovar.pism_name = 'staggered grid coordinate'
    ovar.units = x.units

    # copy time coordinate
    ovar = ods.createVariable('time', t.dtype, ('time'))
    copy_attributes(t, ovar)
    ovar[:] = ts[:]
    ovar.long_name = 'time'
    ovar.units = 's'

    # copy boot bedrock topography
    bvar = bds.variables['topg']
    ovar = ods.createVariable('B', bvar.dtype, ('index1'))
    ovar[:] = bvar[:, :, xcond]
    copy_attributes(bvar, ovar)
    ovar.long_name = 'bed elevation'

    # copy boot ice thickness
    bvar = bds.variables['thk']
    ovar = ods.createVariable('H', bvar.dtype, ('index1'))
    ovar[:] = bvar[:, :, xcond]
    copy_attributes(bvar, ovar)
    ovar.long_name = 'ice thickness'

    # copy effective pressure
    evar = eds.variables['effbwp']
    ovar = ods.createVariable('N', evar.dtype, ('time', 'index1'))
    ovar[:] = evar[tcond, xcond]
    copy_attributes(evar, ovar)
    ovar.long_name = 'effective pressure'

    # copy water sheet thickness
    evar = eds.variables['bwat']
    ovar = ods.createVariable('h', evar.dtype, ('time', 'index1'))
    ovar[:] = evar[tcond, xcond]
    copy_attributes(evar, ovar)
    ovar.long_name = 'water sheet thickness'

    # compute water sheet discharge
    u = eds.variables['bwatvel[0]'][tcond, xcond]
    v = eds.variables['bwatvel[1]'][tcond, xcond]
    ovar = ods.createVariable('q', evar.dtype, ('time', 'index2'))
    ovar[:] = (u**2+v**2)**0.5/(365.0*24*60*60)
    ovar.long_name = 'water sheet discharge'
    ovar.units = 'm^2/s'

    # close datasets
    bds.close()
    eds.close()
    ods.close()


if __name__ == '__main__':
    """Main program, prepare all input files."""

    # create inputs directory if absent
    if not os.path.exists('processed'):
        os.makedirs('processed')

    # postprocess experiments completed
    for exp in exps:
        postprocess(exp)

    # make tar file
    tarname = 'processed/shmip_%s.tar.gz' % auth
    print "Preparing archive %s..." % tarname
    with tarfile.open(tarname, "w:gz") as tar:
        for exp in exps:
            ofilename = 'processed/%s_%s.nc' % (exp.upper(), auth)
            if os.path.isfile(ofilename):
                tar.add(ofilename)

    # copy questionnaire
    shutil.copyfile('questions.rst', 'processed/questions.rst')
