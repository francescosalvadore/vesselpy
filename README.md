# *vesselpy* - Hydrostatic solver

## Installation

[1] Clone the repository:
```
cd $HOME && mkdir VESSELPY && cd VESSELPY
git clone https://github.com/francescosalvadore/vesselpy.git
```

[2] Install requirements (virtual env is recommended):
```
python3 -m venv pyenv
. pyenv/bin/activate
cd vesselpy
pip3 install -r requirements.txt
```

[3] Add vesselpy bin to PATH environment variable (to `.bashrc`):
```
export PATH=$HOME/VESSELPY/vesselpy/bin:$PATH
```

[4] Download and install ParaView:
```
apt install paraview
```

[5] Set ParaView folder in `load_env.sh`, e.g, e.g..:
```
export PARAVIEW_FOLDER="/opt/paraviewopenfoam510/lib/python3.10/site-packages"
```

[6] Run the example:
```
cd .. && mkdir run && cd run
cp -r ../vesselpy/examples/w3/ .
cd w3
vesselpy
```

[7] Check the results:
```
cat output/1/results.dat
cat output/1/stab.dat
```

## Preparing a case

Two input files are needed to start:
* `config.ini` input file
* `stl/<geometry_name>.stl` triangularized (ASCII) geometry

## config.ini

The file is a .ini file including two sections: *general* and *case_\<casename\>* sections.

*general* section includes:
- *gravity*
- atmospheric pressure *p_atmo* (which has no actually no role)
- *refine_factor* which decides how to refine the surface (can be kept as 1 to avoid refinement) 
- *case* which is the name of the case to simulate. The name of the case is used to label the other section to be considered in the run (so that many cases can be kept in the same input file and the *general* section decides which case is run.

An example is given below:
```
[general]
gravity          = 9.81
p_atmo           = 0.   
refine_factor    = 1.
case             = w3
```

*case_\<casename\>* section contains
- *mode*, i.e., the hydrostatic mode which can be on the following:
```
"weight_lcg_tcg", "weight_lcg_heel", "weight_trim_tcg", "weight_trim_heel", "sinkage_trim_heel"
```
corresponding to 3DOF, 2DOF (free trim), 2DOF (free roll), 1DOF, 0DOF solvers.
- *filepath* contains the name of the triangularized geometry in the *stl* folder
- *mass* of the hull
- *rho* water density
- *cog_x*, *cog_y*, *cog_z* center of gravity coordinates
- *cor_x*, *cor_y*, *cor_z* center of (trim and heel) rotation coordinates. Convenient choices usually include CoG or (0,0,0) points.
- *z_water_fixed* height of water plane or *auto* if it is an output
- *alfa_fixed* trim angle or *auto* if it is an output
- *heel_fixed* heel (roll) angle or *auto* if it is an output
- *heel_stability* defines the values of heel large-angle stability heel angles considered in the analysis. The syntax is given by start:stop:step, e.g., 0:80:10 means that the values 0,10,20,30,40,50,60,70,80 will be considered.
- *control_points* defines the list of points which are tracked and for which the final hydrostatic equilibrium positions are reported. The list must be inserted as list of list using square brackets, e.g., *[[0.0, 0.212, 0.0289288],[0.0, 0.0, 0.01]]*

Beware: *mass*, *cog_x*, *cog_y*, *z_water_fixed*, *alfa_fixed*, *heel_fixed* can be inputs, outputs, or optional inputs depending on the selected mode. For example:
- the 3DOF solver requires *mass*, *cog_x*, *cog_y* as inputs but *z_water_fixed*, *alfa_fixed*, *heel_fixed* can be specified as *auto* or giving values which are treated as initial values to speed-up the convergence of their solutions. - - the 0DOF solver requires *alfa_fixed*, *heel_fixed* and *z_water_fixed* as inputs and *mass*, *cog_x*, *cog_y* must be specified as *auto* or are in any case ignored.
- For more details see the paper.

An example is given below:
```
[case_w3]
mode             = weight_lcg_tcg
filepath         = w3_real.stl
mass             = 32.578
rho              = 997.75
cog_x            = 0.697
cog_y            = 0.
cog_z            = 0.143
cor_x            = 0.697
cor_y            = 0.
cor_z            = 0.143
z_water_fixed    = auto
alfa_fixed       = auto
heel_fixed       = auto
heel_stability   = 0:0:10
control_points   = [[0.0, 0.212, 0.0289288],[0.0, 0.0, 0.01]]
```

## stl/<geometry_name>.stl

The triangularized geometry has to provided using ASCII stl format. Triangles must be sufficienctly detailed to describe the geometry.
Moreover, even for planar surfaces, triangles lenght must be limited so that the convergence of hydrostatic algorithm is easier.
Watertightness of the geometry is recommended but not strictly necessary.
