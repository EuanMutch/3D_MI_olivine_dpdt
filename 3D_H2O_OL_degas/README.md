# Codes for 3D diffusion of water from melt inclusion-bearing olivine phenocrysts
3D model code for water loss from olivine hosted melt inclusions with degassing boundary conditions
A code for calculating the analytical solution for the anisotropic diffusive equilibration of a melt inclusion is also included here.


Finite element model codes were run using the legacy version of FEniCS 2019.1.0.

Instructions for installing and downloading can be found here: https://fenicsproject.org/download/archive/

Mesh files were created in GMSH v3.0.0 and were converted to xml files using Dolfin convert. Original .geo files were created using a pygmsh script contained in the folder 'Mesh_gen_scripts'.

Gmsh can be installed and downloaded here: https://gmsh.info/#Download

Here is a description of the included files in the subfolders.

#--------------------------------------------------------------------------------------------
#3D_H2O_OL_degas
3D model code for water loss from olivine hosted melt inclusions with degassing boundary conditions

MI_processes.sh - Wrapper script with input options for the diffusive anisotropy and number of timesteps. 

3D_diffusion_Seguam.py - main diffusion modelling script. 


#===========================================================================

The 3D modelling code can be run in a terminal by typing:

bash MI_processes.sh

Make sure to adjust the file with the appropriate conditions, and select the number of cores you would like to run the code on. 

#=============================================================================


write_h5_mesh.py - script to convert xml mesh files produced by Dolfin into .h5 format for use in parallelised code.

Mesh_gen_scripts - folder with python scripts to generate 3D meshes of olivine

olivine_geom.py - pygmsh script to create olivine faces for testing.

xtal_geometries.py - pygmsh script for creating crystal geometries including volumes and spheres

Sphere_geometry.py - functions for calculating plane sphere interaction.

olivine_mesh_gen_r3_fix.py - script for generating olivine and sphrere meshes. 

Mesh files

OL_mesh1 - lower resolution mesh - 25 micron spatial resolution in olivine crystal and 12.5 micron spatial resolution in and around 30 micron radius melt inclusion.

OL_mesh2 - higher resolution mesh - 15 micron spatial resolution in olivine crystal and 7.5 micron spatial resolution in and around 30 micron radius melt inclusion.


Seguam_degas.csv - degassing curve for Seguam magma calculated using VolatileCalc - 5000 steps in the degassing model. This is used as the external boundary condition in the model.
