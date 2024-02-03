# Script to model 3D diffusion of water from olivine hosted melt inclusions using imported mesh

from dolfin import *
import matplotlib.pyplot as plt
import math as m

import pandas as pd
#import vcalc_cy as vc
import numpy as np
from scipy.interpolate import interp1d
import random
import math as m

import os
import sys
import time
import datetime

#parameters["refinement_algorithm"] = "plaza_with_parent_facets" #refinement
set_log_active(False) #Option to supress outputs

#############################################################################
# PARAMETERS FROM TERMINAL

f_mesh = sys.argv[1] # Mesh file
f_coords = sys.argv[2] # Mesh coordinates information
f_degas = sys.argv[3] # Degassing curve curve file

# set parameters here
f_aniso = float(sys.argv[4]) # Diffusive anisotropy in olivine
f_nts = int(sys.argv[5]) # Number of time steps
f_Kd = float(sys.argv[6])  # Olivine-melt partition coefficient for water
f_dpdt = float(sys.argv[7]) # Magma decompression rate MPa/s
f_h2oi = float(sys.argv[8]) # Initial water content in the melt 
f_h2oe = float(sys.argv[9]) # water content at the end of decompression path based on embayment wt%
T = float(sys.argv[10]) # Temperature (C)

T_K = T + 273.0 # Convert temperature into Kelvin


###############################################################################
#FILES
# Files to be imported and written

# Import mesh here
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), './{0}/{0}.h5'.format(f_mesh), "r")
hdf.read(mesh, "/mesh", False)
domains = MeshFunction('size_t', mesh, mesh.topology().dim())
hdf.read(domains, "/domains")
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")

# Import dataframe of coordinates and concentrations

df_coords = pd.read_csv('./{0}/{1}.csv'.format(f_mesh, f_coords))

X = df_coords['MI_cen_x'].values # For data collection coordinates
Y = df_coords['MI_cen_y'].values
Z = df_coords['MI_cen_z'].values

mi_c = np.array([X[0], Y[0], Z[0]]) # Define centre of melt inclusion based on coordinates
mm = Point(mi_c[0], mi_c[1], mi_c[2])
mm_cell, mm_distance = mesh.bounding_box_tree().compute_closest_entity(mm)

# Define finite element function space
Q = FunctionSpace(mesh, "CG", 1)

###############################################################################
# TIMESTEP AND DEGASSING CALCULATIONS

maxstep = f_nts #100 # maximum number of time steps

# Import degassing curve from file
degas = pd.read_csv(f_degas)

vc_P = degas['P_MPa'].values # Pressure 
vc_h2o = degas['H2O'].values # Water content 

dsteps = len(vc_h2o)

dgci = interp1d(vc_h2o, vc_P)
Ps, Pe = vc_P[0], dgci(f_h2oe)
dp = (Ps - Pe)/maxstep # change in pressure at each step

# Calculate timestep based on decompression rate and change in pressure

dpdt = f_dpdt #MPa/s 0.01 0.16
dt = dp/dpdt
dt_c = Constant(dt) 

# interpolate degassing curve onto new array of pressures if degassing curve != maxstep

if dsteps != maxstep:
    dgc = interp1d(vc_P, vc_h2o)
    Pnew = np.linspace(Ps, Pe, maxstep)
    H2Onew = dgc(Pnew)   # use interpolation function returned by `interp1d`

else:
    H2Onew = vc_h2o

################################################################################
# DIFFUSION COEFFICIENTS

# Olivine diffusion coefficient - H+ with anisotropy - values from Barth et al. (2019)
D_100 = ((9.6e-6)*m.exp(-125000.0/(8.314*T_K)))*(1e12)
D_010 = (1.0/f_aniso)*D_100 # Assume only 10x anisotropy for now. 40 
D_001 = (1.0/f_aniso)*D_100

D0 = as_matrix(((Constant(D_100), Constant(0.0), Constant(0.0)), (Constant(0.0), Constant(D_010), Constant(0.0)), (Constant(0.0), Constant(0.0), Constant(D_001))))

# Calculate diffusivity in the melt inclusion
def Dwtre(Col, Kd):
    Cmelt = (Col/Constant(Kd))/Constant(10000.0)
    Dz = Constant(1.546e-12)*(Cmelt**Constant(3.0)) + Constant(-1.302e-12)*(Cmelt**Constant(2.0)) + Constant(7.442e-12)*(Cmelt) + Constant(1.522e-12)
    return Dz


def D1(C, Kd):
    Dd = Dwtre(C, Kd)*Constant((1e12))
    D1 = as_matrix(((Dd, Constant(0.0), Constant(0.0)), (Constant(0.0), Dd, Constant(0.0)), (Constant(0.0), Constant(0.0), Dd)))
    return D1


#######################################################################################
# MESH TAGGING

cell_dim = mesh.topology().dim()
cmarker = MeshFunction("size_t", mesh, cell_dim, 0)
fmarker = MeshFunction("size_t", mesh, cell_dim - 1, 0)

vmarker = MeshFunction("size_t", mesh, cell_dim - 2, 0)
tmarker = MeshFunction("size_t", mesh, cell_dim - 2, 0)

bmarker = MeshFunction("bool", mesh, cell_dim, False)

# Section here to map tags from gmsh to tags used below

MI_flag = 2
Ol_flag = 4

eb = 3
ib = 1

for c in cells(mesh):
    if domains[c] == MI_flag:
        cmarker[c] = 1
    else:
        cmarker[c] = 0

for f in facets(mesh):
    if boundaries[f] == ib:
        fmarker[f] = 1
    elif boundaries[f] == eb:
        fmarker[f] = 2
    else:
        fmarker[f] = 0


for c in cells(mesh):
    for f in facets(c):
        if fmarker[f] == 2:
            bmarker[c] = True
        elif fmarker[f] == 1:
            bmarker[c] = True

####################################################################################
# CREATE FUNCTION SPACE AND BOUNDARY CONDITIONS

dr_olm = 1.2 # density ratio of olivine/melt

dr_mol = 1.0/1.2 # density ratio of melt/olivine

# Define function spaces
C0 = Function(Q)         
C1 = Function(Q)
S = TestFunction(Q)

# An expression which maps the MeshFunction cmarker to cells
class CellMark(UserExpression):
    def eval_cell(self, values, x, cell):
        values[0] = cmarker[cell.index]
    def value_shape(self):
        return ()
cm = CellMark()

# Create a DG0 Function which has the same values as cmarker
DG0 = FunctionSpace(mesh, "DG", 0)
d_mark = interpolate(cm, DG0)

# Theta method:
# theta = 0.0 is forward Euler
# theta = 0.5 is Crank-Nicholson
# theta = 1.0 is backward Euler
theta = Constant(0.5)
C_mid = theta*C1 + (Constant(1.0)+Constant(-1.0)*theta)*C0

# Set the partition coefficient
Kd = f_Kd  # Olivine/melt partition coefficient from Barth et al., (2019) # 0.001

dofmap = Q.dofmap()
d1_dofs = []
d0_dofs = []

for c in cells(mesh): # compute dofs in the domains
    if cmarker[c] == 1:
        d1_dofs.extend(dofmap.cell_dofs(c.index()))
    else:
        d0_dofs.extend(dofmap.cell_dofs(c.index()))

# unique
d1_dofs = list(set(d1_dofs))
d0_dofs = list(set(d0_dofs))

# Initial conditions for olivine 

# Use degassing curve start point as initial condition

ol_h2o = (H2Onew*10000.0)*Kd
ol_init = ol_h2o[0]
C0.vector()[:] = Constant(ol_init)

# Combined diffusion coefficients for different parts of the mesh
D = D1(C_mid, Kd)*d_mark + D0*(1-d_mark)
rho = ((dr_mol)/Kd)*d_mark + 1.0*(1-d_mark)

F = rho*S*(C1-C0)*dx + dt_c*(dot(rho*D*grad(C_mid), grad(S)))*dx

a, L = lhs(F), rhs(F)
J = derivative(F, C1)

# boundary conditions to start with
Cbc0 = DirichletBC(Q, ol_h2o[0], fmarker, 2)  # fixed boundary conditions on exterior of mesh
Cbcs = [Cbc0]

################################################################################
# TIMESTEPPING

ts = np.empty([0])
mi_cc = np.empty([0])
Pp = np.empty([0])

out_file = File("output_seguam/conc.pvd")

# Time-stepping
i = 0
t = 0
while i < maxstep:
    P = Pnew[i]/1000.0
    ts = np.append(ts, t)
    Pp = np.append(Pp, P*1000.0)

    # changing boundary conditions
    Cbc0 = DirichletBC(Q, Constant(ol_h2o[i]), fmarker, 2)  # fixed boundary conditions on exterior of mesh - changing with degassing curve.
    Cbcs = [Cbc0]
    solve(F==0, C1, bcs=Cbcs, J=J, solver_parameters={'newton_solver':{'linear_solver': 'cg', 'preconditioner': 'hypre_amg'}}) # Solve with CG and algebraic multigrid
    out_file << C1
    C0.assign(C1)

    # Work out C0(mm) for a parallel mesh
    u_eval = C0(*mm) if mm_distance < DOLFIN_EPS else None
    computed_u = mesh.mpi_comm().gather(u_eval, root=0)
    computed_u = mesh.mpi_comm().bcast(computed_u, root=0)
    value = [u for u in computed_u if u is not None][0]
    mi_cc = np.append(mi_cc, value)
    
    if mesh.mpi_comm().rank == 0:
        print('step {0}/{1}'.format(i, maxstep))

    t += dt
    i += 1

print('Model Complete!!!')

####################################################################################
# OUTPUT DATA TO FILES - optional

if mesh.mpi_comm().rank == 0:
    df_mio = pd.DataFrame({'Time': ts, 'P': Pp, 'C_mi': mi_cc/f_Kd, 'H2O_melt': H2Onew, 'H2O_ol_boundary': ol_h2o})

    df_mio['MI_x'] = mi_c[0]
    df_mio['MI_y'] = mi_c[1]
    df_mio['MI_z'] = mi_c[2]

    df_mio['MI_radius'] = df_coords['MI_radius'].values[0]
    df_mio['Ol_length_x'] = df_coords['Ol_x'].values[0]
    df_mio['Ol_length_y'] = df_coords['Ol_y'].values[0]
    df_mio['Ol_length_z'] = df_coords['Ol_z'].values[0]
    df_mio['Tolerance'] = df_coords['tolerance'].values[0]

    df_mio['H2O_initial_{0}'.format(f_aniso)] = f_h2oi
    df_mio['dpdt_{0}'.format(f_aniso)] = dpdt
    df_mio['Time_step_{0}'.format(f_aniso)] = dt
    df_mio['Anisotropy_{0}'.format(f_aniso)] = f_aniso
    df_mio['Kd_{0}'.format(f_aniso)] = f_Kd

    df_mio.to_csv('./{0}_mi_out_aniso_{1}_{2}.csv'.format(f_mesh, f_aniso, f_nts), sep=',')


    # Add parameters to a compiled spreadsheet

    df_coords['C_mi_end_{0}'.format(f_aniso)] = value/f_Kd
    df_coords['Mesh_{0}'.format(f_aniso)] = f_mesh
    df_coords['Date_time_{0}'.format(f_aniso)] = datetime.datetime.now()

    df_coords['H2O_initial_{0}'.format(f_aniso)] = f_h2oi
    df_coords['dpdt_{0}'.format(f_aniso)] = dpdt
    df_coords['Time_step_{0}'.format(f_aniso)] = dt
    df_coords['N_Time_steps_{0}'.format(f_aniso)] = f_nts
    df_coords['Anisotropy_{0}'.format(f_aniso)] = f_aniso
    df_coords['Kd_{0}'.format(f_aniso)] = f_Kd
    df_coords['Max_time_{0}'.format(f_aniso)] = t
    df_coords['P_initial_{0}'.format(f_aniso)] = Pnew[0]
    df_coords['dpdt_actual_{0}'.format(f_aniso)] = Pnew[0]/t

    df_coords.to_csv('./{0}.csv'.format(f_coords), sep=',') 









