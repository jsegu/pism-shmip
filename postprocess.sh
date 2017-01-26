#!/bin/bash

# short author name
auth=jseg
inst="ETH Zürich"
name="Julien Seguinot"

# create directory if missing
mkdir -p $auth

# loop on experiments completed
for exp in a{1..6} b{1..5} e{1..5}
do

    # set title
    echo "Postprocessing experiment ${exp^}..."
    title="PISM experiment ${exp^}."

    # input and output files
    ifile="output/${exp}_extra.nc"
    ofile="$auth/${exp^}_$auth.nc"

    # boot file for topographies
    if [ "${exp:0:1}" == "e" ]
    then
        bfile="input/boot_${exp}.nc"
    else
        bfile="input/boot_sqrt.nc"
    fi

    # select the last time slice in extra file
    ncks --dimension time,-1 --variables effbwp $ifile $ofile --overwrite

    # append boot topographies
    ncks --variables thk,topg $bfile $ofile --append

    # add global attributes
    ncatted -a title,global,c,c,"$title" $ofile --overwrite
    ncatted -a meshtype,global,c,c,"structured" $ofile --overwrite
    ncatted -a dimension,global,c,c,"2D" $ofile --overwrite
    ncatted -a channels_on_edges,global,c,c,"no" $ofile --overwrite
    ncatted -a institution,global,c,c,"$name, $inst" $ofile --overwrite

    # rename variables
    ncrename --variable effbwp,N $ofile --overwrite
    ncrename --variable thk,H $ofile --overwrite
    ncrename --variable topg,B $ofile --overwrite

break
done
