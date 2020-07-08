#!/bin/ksh

#Convert hdf to netcdf, add time variable and rename fill value.

path_data=/data/rd_exchange2/tcandela/STAMM/ressources/VGPM/tmp/
prefix=vgpm.
year0=2019
year1=2019
start=273
end=361


echo "Converting hdf to netcdf"
python $stit/preprocessing/hdf2netcdf.py ${path_data} ${prefix} ${year0} ${year1} ${start} ${end}

echo "Adding time variable"
python $stit/preprocessing/Change_netcdf_time.py ${path_data} ${prefix} ${year0} ${year1} ${start} ${end}

echo "Renaming Fill Value"
python $stit/preprocessing/Rename_FillValue.py ${path_data} ${prefix} ${year0} ${year1} ${start} ${end}
