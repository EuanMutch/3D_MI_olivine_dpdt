import sys
from dolfin import *

f_mesh = sys.argv[1] 

# Import mesh here
mesh = Mesh('./{0}.xml'.format(f_mesh))
domains = MeshFunction('size_t', mesh, './{0}_physical_region.xml'.format(f_mesh))
boundaries = MeshFunction('size_t', mesh, './{0}_facet_region.xml'.format(f_mesh))

hdf = HDF5File(mesh.mpi_comm(), './{0}.h5'.format(f_mesh), "w")
hdf.write(mesh, "/mesh")
hdf.write(domains, "/domains")
hdf.write(boundaries, "/boundaries")
