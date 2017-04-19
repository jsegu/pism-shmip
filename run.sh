#!/bin/bash

# Environment variables
# ---------------------

: ${PISM_DO:=""}                # use "echo" to just display the command
: ${PISM_MPIDO:=""}             # MPI command for parallel run
: ${PISM_EXEC:="pismr"}         # name of PISM executable


# Command-line argument parser
# ----------------------------

# loop on keyword, argument pairs
while [[ $# -gt 0 ]]
do
    case "$1" in
        -e|--exp)
            exp="$2"
            shift
            ;;
        -c|--config)
            config="$2"
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
config=${config:=""}


# Experiment specific set4tings
# ----------------------------

# boot arguments
vert_grid="-Mz 2 -Lz 2000"
case $exp in
    a*|b*)
        boot_args="-i input/boot_sqrt.nc -bootstrap -Mx 403 -My 41 $vert_grid"
        ;;
    c*)
        boot_args="-i output/b5.nc"
        ;;
    e*)
        boot_args="-i input/boot_$exp.nc -bootstrap -Mx 303 -My 56 $vert_grid"
        ;;
    d*|f*)
        echo "Sorry, exp. $exp not implemented yet."
        exit 2
        ;;
esac

# config overrides
case $exp in
    a*|b*|c*|d*)
        conf_file="config_sqrt${config:+_$config}"
        ;;
    e*|f*)
        conf_file="config_valley${config:+_$config}"
        ;;
esac

# basal melt
case $exp in
    a1)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 7.93e-11"
        ;;
    a2)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 1.59e-9"
        ;;
    a3)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 5.79e-9"
        ;;
    a4)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 2.5e-8"
        ;;
    a5)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 4.5e-8"
        ;;
    a6)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 5.79e-7"
        ;;
    a6)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 5.79e-7"
        ;;
    b*)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 7.93e-11"
        melt_args+=" -hydrology_input_to_bed_file input/melt_$exp.nc"
        ;;
    c*)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 7.93e-11"
        melt_args+=" -hydrology_input_to_bed_file input/melt_$exp.nc"
        # PISM does not support non-integer periods, e.g. 1/365 (issue #380).
        #melt_args+=" -hydrology_input_to_bed_period 0.0027397260273972603"
        ;;
    e*)
        melt_args="-hydrology_use_const_bmelt -hydrology_const_bmelt 1.158e-6"
        ;;
esac

# run duration and output
case $exp in
    a*|b*|e*)
        years="5"
        ex_dt=monthly
        ts_dt=daily
        ;;
    c*)
        years="0.082191781"  # 30 days
        ex_dt=hourly
        ts_dt=hourly
        ;;
esac


# Output settings
# ---------------

# run name
mkdir -p output
run="output/$exp${config:+_$config}"

# extra and ts variables
extra_vars=bwat,bwatvel,bwp,bwprel,effbwp,tauc,wallmelt  # diagnostics
#extra_vars+=,hydrobmelt,hydroinput,hydrovelbase_mag  # verification
ts_vars=hydro_ice_free_land_loss,hydro_ice_free_land_loss_cumulative,
ts_vars+=hydro_negative_thickness_gain,hydro_negative_thickness_gain_cumulative


# Actual run
# ----------

# build config file
ncgen $conf_file.cdl -o $conf_file.nc

# run PISM
$PISM_DO $PISM_MPIDO $PISM_EXEC \
    -config_override $conf_file.nc -report_mass_accounting \
    $boot_args $melt_args \
    -y $years -o $run.nc -o_size small \
    -extra_file ${run}_extra.nc -extra_times $ex_dt -extra_vars $extra_vars \
    -ts_file ${run}_ts.nc -ts_times $ts_dt -ts_vars $ts_vars \
    > $run.log 2> $run.err
