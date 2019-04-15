#!/bin/bash

# Parameters
# -------------------------------------------------------------------------------

# PART TO EDIT
# You also have to edit the namelist to specify passive or active and initial positions

# # Passif VGPM with initial positions computed from new forcings
# outdir=/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/run/passif_vgpm/
# namelist=/homelocal/otitaud/gitlab/lmtl-wp4/config/namelists/namelist_glorys2v4_npp_R2018_vgpm_passif

# Active VGPM with initial positions computed from new forcings (mmolC)
# outdir=/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/run/actif_vgpm/
# namelist=/homelocal/otitaud/gitlab/lmtl-wp4/config/namelists/namelist_glorys2v4_npp_R2018_vgpm_actif

# Test
outdir=/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/run/actif_vgpm/test
namelist=/homelocal/otitaud/gitlab/lmtl-wp4/config/namelists/namelist_glorys2v4_npp_R2018_vgpm_actif


## Passif VGPM with reference initial positions
# outdir=/data/FSHML/Tortues/LMTL-WP4/mercatorglorys2v4_1998_2015/run/passif_vgpm_refpos/
# namelist=/homelocal/otitaud/gitlab/lmtl-wp4/config/namelists/namelist_glorys2v4_npp_R2018_vgpm_passif.refpos


# Execution
# -------------------------------------------------------------------------------

init_file=$(grep init_file ${namelist} | awk -F "=" '{print $2}')
echo "Using init file ${init_file}"

mkdir -p ${outdir}

init_file=$(echo ${init_file} | tr -d \')

bname=$(basename ${namelist})
bname=${bname/namelist_/}
outfile=${outdir}/${bname}.nc

cp -v ${namelist} ${outdir}

cp -v ${init_file} ${outdir}

/usr/bin/python2 ../stamm/IBM2D.py $namelist $outfile
res=$?

if [[ "${res}" = "0" ]]; then
    echo "creating figure"
    ../scripts/fig_plot_dispersion.py ${outfile}
else
    echo "STAMM failed"
fi
        
