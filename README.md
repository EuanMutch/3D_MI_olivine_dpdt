# 3D_MI_olivine_dpdt
Model code to support 3D diffusion of water in melt inclusion-bearing olivine phenocrysts

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

water_melt.py - calculating water diffusion in melt using Ni and Zhang (2018).

3D_diffusion_Seguam.py - main diffusion modelling script. 


#===========================================================================

The 3D modelling code can be run in a terminal by typing:

bash MI_processes.sh

Make sure to adjust the file with the appropriate conditions, and select the number of cores you would like to run the code on. 

#=============================================================================

Mesh_gen_scripts - folder with python scripts to generate 3D meshes of olivine

Mesh_gen.sh - bash shell script to generate mesh using olivine_mesh_gen.py and gmsh. Converts the mesh into xml using dolfin. Then converts it into h5 using write_h5_mesh.py

olivine_geom.py - pygmsh script to create olivine faces for testing.

xtal_geometries_2021.py - pygmsh script for creating crystal geometries including volumes and spheres

Sphere_geometry.py - functions for calculating plane sphere interaction.

olivine_mesh_gen.py - script for generating olivine and sphrere meshes. 

olivine_mesh_gen_MC.py - script for generating olivine and sphrere meshes where the melt inclusion size/position and crystal size are randomly selected.

write_h5_mesh.py - script to convert xml mesh files produced by Dolfin into .h5 format for use in parallelised code.

Mesh files

OL_mesh1 - lower resolution mesh - 25 micron spatial resolution in olivine crystal and 12.5 micron spatial resolution in and around 30 micron radius melt inclusion.

Seguam_degas.csv - degassing curve for Seguam magma calculated using VolatileCalc - 5000 steps in the degassing model. This is used as the external boundary condition in the model.

#=============================================================================

Monte_Carlo_scripts - subfolder containing mesh generation scripts in Monte Carlo modelling

olivine_mesh_gen_r5.py - generate 3D olivine mesh with random parametes. Same as olivine_mesh_gen_MC.py. 

1Dp_mesh_gen_r2_asymm_v3.py - generate 1D asymmetric mesh given 3D model input parameters.

1Dp_mesh_gen_r2_symm_v3.py - generate 1D symmetric mesh given 3D model input parameters.

mi_proc_2D_mesh_v11.py - generate 2D asymmetric mesh given 3D model input parameters. Can generate mesh along xy or xz slice.

#================================================================================

 
#----------------------------------------------------------------------------------------------
# Anisotropic_analytical

anisotropic_analytical.py - Python script for calculating the analytical solution for the anisotropic diffusive equilibration of a melt inclusion. This uses equations 14-16 in the main manuscript.

Seguam_degas.csv - degassing curve for Seguam magma calculated using VolatileCalc - 5000 steps in the degassing model. This is used as the external boundary condition in the model.


#----------------------------------------------------------------------------------------------
# V_calc_cyphon

VolatileCalc python and cyphon scripts used to calculate volatile degassing trends.

VolatileCalc.py - VolatileCalc code from Rasmussen et al., (2020).

vcalc_cy.pyx - Cython version of the Volatilecalc. You will need to compile and set this up on your own system. Follow the instructions here for help: https://www.youtube.com/watch?v=Ic1oE6SEOBs&t=244s

setup.py - script to help set up cython version

vcalc_cy.c - c version of cython VolatileCalc code

