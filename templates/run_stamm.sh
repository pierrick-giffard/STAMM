#!/bin/bash

# -------------------------------------------------------------------------------
# Choose directory and namelist
# -------------------------------------------------------------------------------
dir=/home/EXP/GLORYS12/pacifique/active
namelist=${dir}/namelist.txt

# -------------------------------------------------------------------------------
# Execution
# -------------------------------------------------------------------------------
init_file=$(grep init_file ${namelist} | awk -F "=" '{print $2}')
echo "Using init file ${init_file}"


init_file=$(echo ${init_file} | tr -d \')

bname=$(basename ${namelist})
bname=${bname/namelist_/}
outfile=${dir}/${bname}.nc


python /data/rd_exchange2/pgiffard/stamm/src/STAMM.py $namelist $outfile
res=$? 


if [[ "${res}" = "0" ]]; then
    rm -rf out-*
    echo "creating figure"
    python /data/rd_exchange2/pgiffard/stamm/STIT/fig_plot_dispersion.py ${outfile}
else
    echo "STAMM failed"
fi
        
