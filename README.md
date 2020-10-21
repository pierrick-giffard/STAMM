# STAMM

The Sea Turtle Active Movement Model (STAMM) is an Individual Based Model that aims at calculating juvenile sea turtles dispersal. A lagrangian approach is used in which turtles ground velocity is computed as the result of the ocean current velocity and a swimming velocity directed towards favorable habitats (Gaspar and Lalire 2017). Originally based on the lagrangian software ARIANE, its last version uses PARCELS (http://oceanparcels.org/), a set of Python classes and methods for lagrangian simulations. 

Please note that the model is still under development and that you might need to look deeper in the code for specific uses.

# Folder Contents

There are four directories:

- A directory **src** including the source code of STAMM.
- A directory **templates** with examples of namelist, script to run the model, release file,...
- A directory **STIT** containing pre-processing and post-processing tools.
- A directory **LIB** containing libraries used in STIT and in src.

# Set up your python environment

You can create a new environment named _stamm_ and including necessary modules with the following command:

_conda create -n stamm -c conda-forge python=3.6 parcels=2.2.0 spyder basemap basemap-data-hires opencv_

Then, activate this environment with command _conda activate stamm_.

Note: STAMM was tested with Parcels 2.2.0

In your _~/.bashrc_ file, you can add command _export stamm=/data/rd_exchange3/pgiffard/stamm/_ to have stamm directory as shortcut.

# Forcings

## General considerations

STAMM needs 4 variables to run: U and V, respectively the zonal and meridional surface current velocities; T, the ocean surface temperature; and a food proxy (Net Primary Production - NPP for example).

The model accepts regular and non-regular grids (Orca for example), A-grid or C-grid, at any time resolution. PARCELS will take care of spatial and temporal interpolations. The food proxy (NPP for example) can be on a different grid than physics variables (U, V and T). For the moment, STAMM was tested with A-grids only. The food proxy needs to be vertically integrated.

## Download data from CMEMS

You can download data directly on the CMEMS website. We also created a script to help downloading a whole dataset: _STIT/preprocessing/Extract\_CMEMS\_data.py_.

## Files format

Parcels needs data files to satisfy the following conditions:

- They need to have a time variable giving the file date. Files have to be ordered in time. You can use the script _STIT/preprocessing/Change\_netcdf\_time.py_.
- Variables need a FillValue attribute to represent land or cells without data. If Parcels doesn&#39;t find it, please rename it using the script _STIT/preprocessing/Rename\_FillValue.py_ (for example it doesn&#39;t find Hole\_Value).
- PARCELS also needs to read the coordinates of your cells in longitude/latitude: the cells centers for A-grids and f-nodes (upper right corner) for C-grids. It can read them directly in a forcing file or in another mesh file. You can create a mesh\_file following the VGPM example: _STIT/preprocessing/Create\_VGPM\_mesh.py_.

# Release file

Create a release file giving release longitude, latitude and time for each turtle. Follow the example in templates. You can use the script _STIT/create\_init\_pos.py ._ STAMM needs time in days since the 1st of january of the starting year.

# Calibrate model

Calibrate model for your species and for the zone your turtles will move in. To calibrate P0 as the 90th percentile of the NPP distribution in a zone, you can use the script _STIT/preprocessing/P0\_calibration.py_.
You will need to set parameters in namalist as well as in the _src/Species_ folder.

# Namelist

Below you will find a description of all STAMM parameters.

Note: no space should appear between two sections.

**&amp;GENERAL PARAMETERS**

**init\_file** = path to your release file

**nturtles** = number of turtles

**tstep** = computational time step in seconds

**ndays\_simu** = Number of days to simulate

**time\_periodic** = To loop periodically data. It is set to either False, the length of the period in days (since first particle release) or &#39;auto&#39;. In mode &#39;auto&#39;, it will loop automatically after last complete available year. Default: &#39;auto&#39;.

**time\_extrapolation** = To repeat the last datafiles available until the simulation end. Default: False

**t\_output** = Writing time of output file.

**mode** = Simulation mode: &#39;passive&#39; (ocean currents only) or &#39;active&#39; (ocean currents plus swimming velocity). Default: active.

**key\_alltracers** = boolean, if False doesn&#39;t use T and NPP for diagnostics. Default: True.

**adv\_scheme** = advection scheme: &#39;RK4&#39; or &#39;Euler&#39;. Default: RK4

**ystart** = integer, simulation starting year.

/

**&amp;ACTIVE MODE**

**species** = &#39;Leatherback&#39; or &#39;Loggerhead&#39;.

**alpha** = Scaling parameter for the Von Mises distribution in m-1.

**P0** = Scaling factor used to calculate food habitat expressed in the same units as NPP.

**grad\_dx** = distance from turtle (in meters) to which habitat gradients are computed.

**cold\_death** = boolean, if True turtles are deleted after cold\_resistance days under Tmin. Default: False

**cold\_resistance** = integer, number of days turtles can resist under Tmin. Default: 30.

**growth** = &#39;VGBF&#39; or &#39;Gompertz&#39;. To use a von Bertalanffy function for growth or the Gompertz model with habitat dependance (only for Leatherbacks).

**SCL0** = &#39;float&#39;. Turtles initial SCL.

**tactic\_factor** = float [0, 1]. Proportion of current swimming angle taken into account in final swimming direction. 1 for no memory effect. Default: 1.

**frenzy** = boolean. To simulate frenzy swimming. Turtles will swim during several days towards the same direction (define parameters in src/Species). Default: False

**wave\_swim** = boolean. So that turtles swim against wave at the beginning of the simulation. Wave direction is calculated based on Stokes drift. Default: False

**dt\_swim** = float. Swimming velocity is considered constant over dt\_swim seconds. Default: 86400

**/**

**&amp;GRIDPARAM**

**periodicBC** = boolean, set to True for East/West periodicity. Default: True.

**halo** = boolean. Set to True if your data already has a halo (for example Orca grids already has a repetition of last/first columns). Default: False.

/

**&amp;FORCINGS FILES**

**U\_dir** = data directory of the zonal velocity

**U\_suffix** = files suffix

**V\_dir** = data directory of the meridional velocity

**V\_suffix** = files suffix

**T\_dir** = data directory of the temperature

**T\_suffix** = files suffix

**food\_dir** = data directory of the food proxy.

**food\_suffix** = files suffix

/

**&amp;PHYSICS**

**mesh\_phy** = File from which lon/lat coordinates will be extracted for physics (U, V and T). Can be one of the data files for A-grids. For a C-grid you need to provide f-nodes (cells corners).

**grid\_phy** = &#39;C&#39; for C-grid, &#39;A&#39; for A-grid (centered grid). The interpolation scheme is nearest neighbour for C-grid tracers.

**lon\_phy** = name of longitude variable in mesh\_phy (f-node if C-grid)

**lat\_phy** = name of latitude variable in mesh\_phy (f-node if C-grid)

**lon\_T** = Only for C-grid. Name of temperature longitude variable in mesh\_phy

**lat\_T** = Only for C-grid. Name of temperature latitude variable in mesh\_phy

**time\_var\_phy** = name of time variable in your physical files

**U\_var** = name of zonal velocity variable

**V\_var** = name of meridional velocity variable

**T\_var** = name of temperature variable

/

**&amp;FOOD PROXY**

**mesh\_food** = File from which lon/lat coordinates will be extracted for food proxy. Can be one of the data files for A-grids.

**lon\_food** = name of longitude variable in mesh\_food

**lat\_food** = name of latitude variable in mesh\_food

**time\_var\_food** = name of time variable in your food proxy files

**food\_var** = name of food proxy variable

/

**&amp;WAVES**

**wave\_dir** = data directory for waves (Stokes drift) if wave_swim is True.

**wave\_suffix** = files suffix

**Ust\_var** = name of zonal Stokes drift

**Vst\_var** = name of meridional Stokes drift

**/**

# Execute STAMM

Execute command _python src/STAMM.py namelist outfile.nc_. STAMM will create a netcdf file with turtles trajectories and with many other variables such as habitat, temperature, SCL,...
You can also use script _src/run\_stamm.sh_ to run STAMM, and then to execute some post-processing diagnosis. 


# Analyse and Visualize results with STIT

You will find a few scripts to visualize the outputs in _STIT/postprocessing_.

# General remarks

## Forcings Selection

STAMM is sensitive to forcings selection. In the namelist, you will provide a directory containing one or several variables. For each variable, a suffix is also given. Take care that one suffix matches only with the desired variable. You can also create symbolic links pointing to your files.

If it doesn&#39;t work for your files, try to modify your file names or the function forcing\_list in Funtions.py.



## Land mask

In STAMM, all variables are set to 0 on land so that habitat gradient is directed towards ocean and turtles avoid beaching. Land corresponds to the physical land mask, and more precisely a turtle is on land when the 4 surrounding points are masked points.
Physical land mask can be different from food land mask, but it can generate problems. Most turtles won&#39;t enter into the food land mask because habitat is set to 0. If they do so, they could be &#39;lost&#39; because there won&#39;t be any habitat gradient because habitat will be 0 everywhere.

## Escape beaching
A kernel was created so that turtle do not beach. If a turtle is about to beach,  it stays at the same position. At next time step, the parameter grad\_dx is set at the model grid resolution to avoid calculating ggradient further than 1 cell. In most cases the turtle will swim towards ocean. If it doesn't, at the next timestep it will take the opposite direction.

## Disabled turtles

Turtles can beach in case of convergence or in case the hazard on swimming direction lead them towards land. This shouldn&#39;t happen too often since variables are set to 0 on land, and hence habitat is also 0 on land. The habitat gradient is then directed towards the open ocean. In case turtles beach, they are sent back to their last position. They are disabled (don&#39;t move anymore) if:

- they pass through the northern boundary of domain.
- periodicBC is set to False and they reach the east-west boundary.
- they beach more than 50 times in a row.
- cold\_mortality is set to True and they stay too long in a cold area.

The variable active can help you to visualize disabled turtles.

## Time step

The choice of the simulation time step depends on the advection scheme used, on velocities experimented by turtles and on grid resolution. For Euler advection scheme, the Courant–Friedrichs–Lewy condition has to be respected in the whole area: u \* dt / dx has to be less than 1 in the whole domain. u is the maximum velocity turtles will experiment, dt is the simulation time step and dx is the spatial grid resolution.

The Runge-Kutta scheme is more stable, but has no criteria as simple as for Euler scheme. Hence, and for security reasons, it is possible to use the Courant criteria also for Runge-Kutta 4.

## Computational time and Memory

The computational time is not proportional to time step. Hence, a great part of computational time is dedicated to access data from disk (operation done each day for daily datasets). Also, it is quite fast to calculate turtles positions evolution since it is done in C (up to 10 000 turtles). Indeed, Parcels translates all python kernels to C language in order to speed up execution.

Field chunking is available for large datasets, to load data that would never fit into memory in one go. In STAMM chunksize is set to False for 1/4° resolution 2D fields, to &#39;auto&#39; for 2D fields of higher resolution and to a fixed value (256) for 3D fields. This fixed value might to be optimal in all situations.

Note that using 3D data is a way slower than using surface data only.

## tstep and t\_output

If the writing time step (t\_output) is smaller than the computational time step (tstep), the actual computational time step will be t\_output. On the contrary, if tstep is smaller than t\_output and if they aren&#39;t multiple, Parcels will define a smaller time step multiple of tstep and t\_output.

Also, it&#39;s worth mention that outputs aren&#39;t means over multiple time steps but are only the value at a certain time. For example if tstep =1h and t\_output = 24h, then output variables won&#39;t be the daily mean but only the value each 24h.

# Notes for developers

## Adding a new sea turtle specie

Currently, STAMM is only prepared for Leatherback and Loggerhead sea turtles. To add a new specie, you will need:

1. Create a new parameter file in _src/Species_ using the other ones as templates.
2. In _src/Functions.py_, function _Initialisation_, add an _elif_ for your species to select your parameter file. You also need to import it at the beginning of the file.
3. In namelist, choose your species.

## Adding a new argument to the namelist

You can add argument in the namelist. You only need to add it to the dictionary in _LIB/IOlib.py_

## Adding behaviours to turtles

It is quite easy to add behaviours to turtles. To do so, you can add a function to _Active\_kernels.py_, _Passive\_kernels.p_y or _Advection\_kernel.py._ However, parcels will pass your functions to C code. For this reason, your functions have to use only basic Python: no for loops, no extern function call, …

Then, in Functions.py you need to add your function to the kernel\_list so that it become in kernel used during execution. At each time step, all kernels are applied to each particle one by one. That&#39;s why the order is very important.

In kernels, you can print text as in python, but to print particle attribute you need to use _print('%f'%particle.lon)_.

## C-grid
For now, STAMM was tested with A-grids (centered grids). It should work with C-grids although it hasn't been tested yet.


## Initial tracers sampling
At t=0, all tracers are set to 0. Indeed, defining temperature as _T = Variable('T', initial=fieldset.T)_ is too long to initialize because it is executed in Scipy mode. An other possibility is to execute a sampling kernel for one timestep: _pset.execute(SampleTracers, dt=0)_ but it works only if all particles are released at the same time because of deferred_load (only 3 time steps are loaded at the same time).

