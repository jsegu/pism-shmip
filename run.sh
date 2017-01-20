#!/bin/bash

# executables
: ${PISM_DO:=""}                # use "echo" to just display the command
: ${PISM_MPIDO:="mpiexec -n 8"} # MPI command for parallel run
: ${PISM_EXEC:="pismr"}         # name of PISM executable

# model parameters
: ${RUN:="run"}                 # output files prefix
: ${Y:=1}                       # duration in years

# build config file
ncgen config.cdl -o config.nc

# extra variables
extra_vars=bwat,bwatvel,bwp,bwprel,effbwp,wallmelt  # diagnostics
#extra_vars+=,hydrobmelt,hydroinput,hydrovelbase_mag  # verification

# run PISM
$PISM_DO $PISM_MPIDO $PISM_EXEC \
    -config_override config.nc -report_mass_accounting \
    -i boot_sqrt.nc -bootstrap -o $RUN.nc -o_size small \
    -Mx 403 -My 41 -Mz 2 -Lz 2000 -y $Y \
    -extra_file $RUN-extra.nc -extra_times daily -extra_vars $extra_vars \
    > $RUN.log 2> $RUN.err &
