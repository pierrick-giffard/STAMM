#!/bin/bash

# -------------------------------------------------------------------------------
# Choose directory and namelist
# -------------------------------------------------------------------------------
outdir=./outputs
namelist=./namelist

animation=True

# -------------------------------------------------------------------------------
# Execution
# -------------------------------------------------------------------------------
rm -rf ${outdir}/out-*

init_file=$(grep init_file ${namelist} | awk -F "=" '{print $2}')
init_file=$(echo ${init_file} | tr -d \')

mkdir -p ${outdir}

bname=$(basename ${namelist})
bname=${bname/namelist_/}
outfile=${outdir}/${bname}.nc


python $stamm/src/STAMM.py $namelist $outfile
res=$? 


if [[ "${res}" = "0" ]]; then
    echo "Computing statistics"
    python $stit/postprocessing/Quick_stats.py ${outfile} $namelist    

    echo "Creating figure"
    python $stit/postprocessing/fig_plot_dispersion.py ${outfile}

    if [[ ${animation} = True ]]; then
	echo "Creating animations"
	python $stit/postprocessing/plot_animation.py ${outfile} $namelist
    fi
else
    echo "STAMM failed"
fi
        
