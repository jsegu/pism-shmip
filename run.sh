#!/bin/bash

# executables
: ${PISM_DO:=""}                # use "echo" to just display the command
: ${PISM_MPIDO:=""}             # MPI command for parallel run
: ${PISM_EXEC:="pismr"}         # name of PISM executable

# loop on keyword, argument pairs
while [[ $# -gt 0 ]]
do
    case "$1" in
        -e|--exp)
            exp="$2"
            shift
            ;;
        -y|--years)
            years="$2"
            shift
            ;;
        *)
            echo "Unknown option $1. Exiting."
            exit 0
            ;;
    esac
    shift
done

# default arguments
exp=${exp:="a1"}
years=${years:="5"}

# run name
mkdir -p output
run="output/$exp"

# extra and ts variables
extra_vars=bwat,bwatvel,bwp,bwprel,effbwp,wallmelt  # diagnostics
#extra_vars+=,hydrobmelt,hydroinput,hydrovelbase_mag  # verification
ts_vars=hydro_ice_free_land_loss,hydro_ice_free_land_loss_cumulative,
ts_vars+=hydro_negative_thickness_gain,hydro_negative_thickness_gain_cumulative

# build config file
ncgen config.cdl -o config.nc

# run PISM
$PISM_DO $PISM_MPIDO $PISM_EXEC \
    -config_override config.nc -report_mass_accounting \
    -i input/boot_sqrt.nc -bootstrap -hydrology_bmelt_file input/melt_$exp.nc \
    -Mx 403 -My 41 -Mz 2 -Lz 2000 -y $years -o $run.nc -o_size small \
    -extra_file ${run}_extra.nc -extra_times monthly -extra_vars $extra_vars \
    -ts_file ${run}_ts.nc -ts_times daily -ts_vars $ts_vars \
    > $run.log 2> $run.err
