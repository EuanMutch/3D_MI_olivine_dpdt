#/bin/bash

# wrapper script

file_name1="OL_mesh1" # Name of mesh files
file_name2="OL_MI_info_1" # Name of mesh information file
file_degas="Seguam_degas.csv"


aniso=15.6 #diffusion anisotropy
nts=100 # Number of timesteps
kd=0.000459 #Olivine-melt partition coefficient for water
dpdt=0.02 # Magma decompression rate MPa/s
h2oi=4.2 # Initial water content in the melt
h2oe=0.3 # Water content at the end of decompression path based on embayment wt%
temp=1070.0 #Temperature (C)

mpirun -np 18 python3 3D_diffusion_Seguam.py ${file_name1} ${file_name2} ${file_degas} ${aniso} ${nts} ${kd} ${dpdt} ${h2oi} ${h2oe} ${temp} # Currently being run on 18 cores


