# HiPSTAR_overset_meshing_tool
Meshing tool to generate overset meshes

Prerequisites: 

 - Add HiPSTAR PLATUS tool to python path

```sh
export PYTHONPATH="absolutePathToHipstar/hipstar/PLATUS3/fortran/lib/"
export PYTHONPATH="absolutePathToHipstar/hipstar/PLATUS3/:$PYTHONPATH"
export PYTHONPATH="absolutePathToHipstar/hipstar/PLATUS3/meshing_tools/"
```
 - python 3.5 or newer versions


# How to generate a mesh?

 1. Run 01-surface_distribution.py
 2. Run 02a-mesh\_prepare.py or 02b-mesh\_prepare\_fishtail.py, according to the type of mesh
  - O-grid with one or multiple background blocks
  - O-grid with near-wake refinemnent and background block


# FILES

## 01-surface_distribution.py
Generate a point distribution on the airfoil/blade.

INPUT PARAMETERS:
 - profile : name of airfoil data in geom folder (airfoilname)
 - n : number of points around the airfoil
 - smooth_iter : parameter smoothing 2nd derivatives (in the range 1e4 - 1e6)
 - chord : chord length
 - TEorigin : (boolean) locate the origin of coordinates at the trailing edge or the default location

INPUT FILE:
 - geom/airfoilname.dat

OUTPUT FILE:
 - geom_fine/arifoilname.dat


## 02a-mesh_prepare.py
Generate an o-grid mesh overlappling cartesian background meshes.

Use of output of 01-surface_distribution.py refined airfoil geometry.

INPUT PARAMETERS:
 - foil_path : name of airfoil data in geom\_fine folder (airfoilname)
 - case_name : name of the mesh folder
 - AoA : angle of attack to rotate the airfoil (deg)
 - bckg_n : number of desired background blocks
 - dy_rel : first cell size in normal direction (size of first cell)
 - bl_exp : o-grid expansion ratio in normal direction
 - ogrid_ny : number of points in the normal directio around o-grid
 - nk : number of points in the spanwise direction (to estimate the 3D mesh size)
 - In case of better controlling the expansion ratio of the boundary layer o-grid mesh:
  - bl\_exp\_i : array/list with sequential spansion ratio
  - ogrid\_ny\_i : array/list with sequential number of layers
 - offset : vertical offset in case of desiring moving the airfoil vertically
 - ds\_mesh\_rel : background mesh cell size
 - blks\_in, blks\_out, blks\_bottom, blks\_top : inlet, outlet, bottom and top distance to origin of coordinates
 - In case of requiring extra stretching towards bigger domain at the inlet, outlet, top and/or bottom:
  - stretch\_inlet, stretch\_outlet, stretch\_top, stretch\_bottom : Boolean
  - stretch\_len\_inlet, stretch\_len\_outlet, stretch\_len\_top, stretch\_len\_bottom : extra stretching distance
  - expand\_ratio\_inlet, expand\_ratio\_outlet, expand\_ratio\_top, expand\_ratio\_bottom : expansion ratio
 - In case of more than one background block:
  - blks\_in\_i : array/list with x positions for extra blocks
  - blks\_out\_i : array/list with x positions for extra blocks
  - blks\_bottom\_i : array/list with y positions for extra blocks
  - blks\_top\_i : array/list with y positions for extra blocks

INPUT FILE:
 - geom_fine/arifoilname.dat

OUTPUT FILE:
 - meshes/meshname/camb.dat
 - meshes/meshname/z_r_grid_x.dat
 - meshes/meshname/meshname_info.txt
 - meshes/meshname/mesh.png

## 02b-mesh_prepare_fishtail.py
Generate an o-grid mesh overlappling cartesian background and a near-wake "fishtail" refinement region.

Use of output of 01-surface_distribution.py refined airfoil geometry.

INPUT PARAMETERS:
 - foil_path : name of airfoil data in geom\_fine folder (airfoilname)
 - case_name : name of the mesh folder
 - AoA : angle of attack to rotate the airfoil (deg)
 - bckg_n : number of desired background blocks
 - dy_rel : first cell size in normal direction (size of first cell)
 - bl_exp : o-grid expansion ratio in normal direction
 - ogrid_ny : number of points in the normal directio around o-grid
 - nk : number of points in the spanwise direction (to estimate the 3D mesh size)
 - In case of better controlling the expansion ratio of the boundary layer o-grid mesh:
  - bl\_exp\_i : array/list with sequential spansion ratio
  - ogrid\_ny\_i : array/list with sequential number of layers
 - offset : vertical offset in case of desiring moving the airfoil vertically
 - wake\_stretch : wake stretching law
 - x\_wake : x coordinate of last wake cell in the horizontal direction
 - slope\_wake\_press : limit angle to create the wake, pressure side
 - slope\_wake\_suct : limit angle to create the wake, suction side
 - ds\_max\_wake : maximum size of wake cell
 - ds\_mesh\_rel : background mesh cell size
 - blks\_in, blks\_out, blks\_bottom, blks\_top : inlet, outlet, bottom and top distance to origin of coordinates
 - In case of requiring extra stretching towards bigger domain at the inlet, outlet, top and/or bottom:
  - stretch\_inlet, stretch\_outlet, stretch\_top, stretch\_bottom : Boolean
  - stretch\_len\_inlet, stretch\_len\_outlet, stretch\_len\_top, stretch\_len\_bottom : extra stretching distance
  - expand\_ratio\_inlet, expand\_ratio\_outlet, expand\_ratio\_top, expand\_ratio\_bottom : expansion ratio

INPUT FILE:
 - geom_fine/arifoilname.dat

OUTPUT FILE:
 - meshes/meshname/camb.dat
 - meshes/meshname/z_r_grid_x.dat
 - meshes/meshname/meshname_info.txt
 - meshes/meshname/mesh.png

## 04-mesh_plot.py

Visualization routine for the generated mesh.

INPUT PARAMETERS: included in the command line when running the routine
 - stp : step to mesh plot every stp lines
 - meshname : name of the mesh (folder inside meshes directory, i.e. SAMPLE)
 - monitorsfile : name of monitors file to be included (i.e. Monitors_sample.dat)

EXAMPLE:
```python
python3 04-mesh_plot.py --meshname=SAMPLE --monitorsfile=Monitors\_sample.dat --stp=5
```
Shown SAMPLE mesh, with the corresponding monitor points. The mesh is displayed every 5 grid lines.

# INPUT FILES STRUCTURE

## geom/airfoilname.dat

 - Two columns (x and y coordinates) separated by spaces.
 - No header in the columns.
 - Clockwise turn around the airfoil surface starting from the trailing edge

## Monitors.dat

 - Three columns (x, y, and z coordinates) separated by spaces
 - Header in the columns: x y z

# FOLDERS

## geom

Files extension .dat where the variables x and y are saved in two columns separated by a space.

 - The points around the airfoil go from trailing edge trailing towards the leading edge passing through the pressure side, and then from the leading edge towards the trailing edge passing through the suction side. 
 - Airfoil chord normalized to 1.

## geom_fine

Refined geometry as output of execution of 01-surface_distribution.py

## meshes

Mesh folders with the output of the mesh generation . Each folder contains:

 - camb.dat : camber datapoints
 - z_r_grid_x.dat : mesh grid data file
 - xxx.txt : text file with information on the mesh properties
 - mesh.png : image of the resulting mesh
