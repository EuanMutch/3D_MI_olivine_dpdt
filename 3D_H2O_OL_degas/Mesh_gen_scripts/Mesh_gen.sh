#/bin/bash

# wrapper script for generating mesh and then testing position of melt inclusion


file_name1="OL_mesh3" # Mesh name
file_name2="OL_MI_info_3" # Accompanying file with mesh information


# Define parameters for mesh

ol_z=500.0 # Length of olivine z axis
mi_r=30.0 # Melt inclusion radius
mi_x=0.0 # melt inclusion location along [100]
mi_y=0.0 # melt inclusion location along [010]
mi_z=0.0 # melt inclusion location along [001]

mesh_lres=12.5 # mesh resolution in olivine and farfield
mesh_hres=3.0 # mesh resolution near melt inclusion

python3 olivine_mesh_gen.py $file_name1 $file_name2 $mesh_lres $mesh_hres $ol_z $mi_r $mi_x $mi_y $mi_z

# If you want to have melt inclusion and crystal size randomly generated use:

#python3 olivine_mesh_gen_MC.py $file_name1 $file_name2

# Use gmsh to create mesh

gmsh -3 "${file_name1}.geo" #2>&1 >/dev/null ; then

echo 'Mesh Generation Successful!'

# convert mesh for dolfin 

dolfin-convert "${file_name1}.msh" "${file_name1}.xml"

# Write mesh as h5 file

python3 write_h5_mesh.py ${file_name1}












