#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 16:48:41 2019

@author: ejfm2
"""

# Script to generate olivine mesh with melt inclusion 

# Generate random number for position of melt inclusion in olivine 

# Need to create coordinate system for olivine - consider only doing in northern hemisphere of crystal

import pygmsh
import numpy as np
import meshio
import pandas as pd
import random
import sys

import xtal_geometries_2021 as xtal
import olivine_geom as olgeom
import Sphere_geometry as sphg

f_mesh = sys.argv[1]
#f_coords = sys.argv[2]
f_dat = sys.argv[2]
mesh_name = sys.argv[3]
#dat_filt = sys.argv[4]
f_axis1 = sys.argv[4]
f_axis2 = sys.argv[5]

df_dat = pd.read_csv('{0}.csv'.format(f_dat))

# filter dataframe based on mesh name

df_f = df_dat[df_dat['Mesh'] == 'OL_mesh_{0}'.format(mesh_name)]

ol_z = df_f['Ol_z'].values[0]
ol_x = df_f['Ol_x'].values[0]
ol_y = df_f['Ol_x'].values[0]

ol_zs = 250.0 # z axis in initial construction

sf = ol_z/ol_zs

mi_x = df_f['MI_cen_x'].values[0]
mi_y = df_f['MI_cen_y'].values[0]
mi_z = df_f['MI_cen_z'].values[0]

mi_r = df_f['MI_radius'].values[0]

h2oi = df_f['H2O_initial'].values[0]
#h2oe = df_f['C_mi_end'].values[0]

tol = mi_r/10.0  # edge tolerance
mesh_lres = ol_z/60.0
mesh_hres = mi_r/12.0

# Add conditions for the maximum number of edges 

if f_axis1 == "x" and f_axis2 == "z":
    if mi_y > 0 and mi_y < 15.4053131973*sf:
        n_edges = 8
    if  mi_y > 15.4053131973*sf and mi_y < 46.3202762831*sf:
        n_edges = 8
    if  mi_y > 46.3202762831*sf and mi_y < 82.5153131984*sf:
        n_edges = 4
    if  mi_y > 82.5153131984*sf and mi_y < ol_y:
        n_edges = 4


if f_axis1 == "x" and f_axis2 == "y":
    if mi_z > 0 and mi_z < 121.445412323*sf:
        n_edges = 10
    if  mi_z > 121.445412323*sf and mi_z < 171.279952571*sf:
        n_edges = 10
    if  mi_z > 171.279952571*sf and mi_z < 186.569889325*sf:
        n_edges = 6
    if  mi_z > 186.569889325*sf and mi_z < 213.73673014*sf:
        n_edges = 8
    if  mi_z > 213.73673014*sf and mi_z < ol_z:
        n_edges = 4


# Recalculate melt inclusion position for rotation

if f_axis1 == "x" and f_axis2 == "z":

    mi_y = df_f['MI_cen_z'].values[0]*-1.0
    mi_z = df_f['MI_cen_y'].values[0]    

elif f_axis1 == "z" and f_axis2 == "x": 

    mi_y = df_f['MI_cen_z'].values[0]*-1.0
    mi_z = df_f['MI_cen_y'].values[0]

# Create plane through melt inclusion

# generate sectioning plane in xy

mi_pc = np.array([mi_x, mi_y, mi_z])

mi_px_pve = np.array([mi_x + 4*ol_z, mi_y, mi_z])
mi_pz_pve = np.array([mi_x, mi_y, mi_z + 4*ol_z])
mi_py_pve = np.array([mi_x, mi_y + 4*ol_z, mi_z])

mi_px_nve = np.array([mi_x - 4*ol_z, mi_y, mi_z])
mi_pz_nve = np.array([mi_x, mi_y, mi_z - 4*ol_z])
mi_py_nve = np.array([mi_x, mi_y - 4*ol_z, mi_z])

mi_px_pve2 = np.array([mi_x + mi_r, mi_y, mi_z])
mi_pz_pve2 = np.array([mi_x, mi_y, mi_z + mi_r])
mi_py_pve2 = np.array([mi_x, mi_y + mi_r, mi_z])

mi_px_nve2 = np.array([mi_x - mi_r, mi_y, mi_z])
mi_pz_nve2 = np.array([mi_x, mi_y, mi_z - mi_r])
mi_py_nve2 = np.array([mi_x, mi_y - mi_r, mi_z])

# xy plane

mi_p1r_pve = mi_px_pve
mi_p1r_nve = mi_px_nve

mi_p1r_pve2 = mi_px_pve2
mi_p1r_nve2 = mi_px_nve2


mi_p2r_pve = mi_py_pve
mi_p2r_nve = mi_py_nve

mi_p2r_pve2 = mi_py_pve2
mi_p2r_nve2 = mi_py_nve2

#if f_axis1 == "x":
#    mi_p1r_pve = mi_px_pve
#    mi_p1r_nve = mi_px_nve

#    mi_p1r_pve2 = mi_px_pve2
#    mi_p1r_nve2 = mi_px_nve2

#elif f_axis1 =="y":
#    mi_p1r_pve = mi_py_pve
#    mi_p1r_nve = mi_py_nve

#    mi_p1r_pve2 = mi_py_pve2
#    mi_p1r_nve2 = mi_py_nve2

#elif f_axis1 =="z":
#    mi_p1r_pve = mi_pz_pve
#    mi_p1r_nve = mi_pz_nve

#    mi_p1r_pve2 = mi_pz_pve2
#    mi_p1r_nve2 = mi_pz_nve2

#if f_axis2 == "x":
#    mi_p2r_pve = mi_px_pve
#    mi_p2r_nve = mi_px_nve

#    mi_p2r_pve2 = mi_px_pve2
#    mi_p2r_nve2 = mi_px_nve2

#elif f_axis2 =="y":
#    mi_p2r_pve = mi_py_pve
#    mi_p2r_nve = mi_py_nve

#    mi_p2r_pve2 = mi_py_pve2
#    mi_p2r_nve2 = mi_py_nve2

#elif f_axis2 =="z":
#    mi_p2r_pve = mi_pz_pve
#    mi_p2r_nve = mi_pz_nve

#    mi_p2r_pve2 = mi_pz_pve2
#    mi_p2r_nve2 = mi_pz_nve2

geom = pygmsh.built_in.Geometry()

geom.add_raw_code('SetFactory("OpenCASCADE");')

# Need to generate coordinate reference frame to position melt inclusion

# define olivine size - scale along z axis (positive direction) - in microns 

# generate central point o melt inclusion

# Test location of melt inclusion 
# Intersection with crystal faces and in interior of crystal

# Run tests

mi_centre = [mi_x, mi_y, mi_z]

#mi1 = xtal.sphere2(geom, mi_centre, mi_r, mi_r, mi_r, mesh_hres)

#mi_s = mi1.surface()

#mi_sl = geom.add_surface_loop(mi_s)
#mi_ps = geom.add_physical(mi_s)

#mi_v = geom.add_volume(mi_sl)
#mi_pv = geom.add_physical(mi_v)


geom.add_raw_code('Sphere(1) = {{{0}, {1}, {2}, {3}, -Pi/2, Pi/2, 2*Pi}};'.format(mi_x, mi_y, mi_z, mi_r))


# Create melt inclusion circle first

#p_c0 = geom.add_point([mi_x, mi_y, mi_z], mesh_hres)

#p_c1 = geom.add_point(mi_p1r_pve2, mesh_hres)
#p_c2 = geom.add_point(mi_p2r_pve2, mesh_hres)
#p_c3 = geom.add_point(mi_p1r_nve2, mesh_hres)
#p_c4 = geom.add_point(mi_p2r_nve2, mesh_hres)

#mi_c1 = geom.add_circle_arc(p_c1, p_c0, p_c2)
#mi_c2 = geom.add_circle_arc(p_c2, p_c0, p_c3)
#mi_c3 = geom.add_circle_arc(p_c3, p_c0, p_c4)
#mi_c4 = geom.add_circle_arc(p_c4, p_c0, p_c1)

#ll1 = [mi_c1, mi_c2, mi_c3, mi_c4]

#mi_ll1 = geom.add_line_loop(ll1)

#mi_s1 = geom.add_surface(mi_ll1)

#ol_ps = geom.add_physical(ll1)
#ol_ps = geom.add_physical(mi_s1)

# Create olivine mesh

ol1 = xtal.olivine(geom, ol_zs, mesh_lres)
ol_s = ol1.surfaces()

ol_sl = geom.add_surface_loop(ol_s)
ol_ps = geom.add_physical(ol_s)

#ol_v = geom.add_volume(ol_sl, holes= [mi_sl])
ol_v = geom.add_volume(ol_sl)
ol_pv = geom.add_physical(ol_v)

# scale olivine mesh

if f_axis1 == "x" and f_axis2 == "z":
    geom.add_raw_code('Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} {Volume{99}; }')
elif f_axis1 == "z" and f_axis2 == "x":
    geom.add_raw_code('Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} {Volume{99}; }')

geom.add_raw_code('Dilate {{{{0, 0, 0}}, {{{0}, {0}, {0}}}}} {{Volume{{99}}; }}'.format(sf)) 

# Create sectioning surface - create plane through centre of melt inclusion

sp1 = geom.add_point(mi_p1r_pve, mesh_hres)
sp2 = geom.add_point(mi_p2r_pve, mesh_hres)
sp3 = geom.add_point(mi_p1r_nve, mesh_hres)
sp4 = geom.add_point(mi_p2r_nve, mesh_hres)

# Add lines

sp_l1 = geom.add_line(sp1, sp2)
sp_l2 = geom.add_line(sp2, sp3)
sp_l3 = geom.add_line(sp3, sp4)
sp_l4 = geom.add_line(sp4, sp1)

sp_ll1 = geom.add_line_loop([sp_l1, sp_l2, sp_l3, sp_l4])

sp_s1 = geom.add_surface(sp_ll1)

# Conduct Boolean Intersection 

geom.add_raw_code('BooleanIntersection{ Volume{99}; Volume{1}; Delete; }{ Surface{105}; Delete; }')

# Add if statement depending on whether xy or xz

# Add new line loops

geom.add_raw_code('Line Loop(108) = {{1:{0}}};'.format(n_edges))
geom.add_raw_code('Line Loop(109) = {{{0}:{1}}};'.format(n_edges + 1, n_edges + 2))

geom.add_raw_code('Plane Surface(3) = {108, 109};')
geom.add_raw_code('Physical Line(3) = {{1:{0}}};'.format(n_edges))
geom.add_raw_code('Physical Line(4) = {{{0}:{1}}};'.format(n_edges + 1, n_edges + 2))
geom.add_raw_code('Physical Surface(5) = {1};')
geom.add_raw_code('Physical Surface(6) = {2};')


geom.add_raw_code('Field[1] = Attractor;')
geom.add_raw_code('Field[1].EdgesList = {{{0}:{1}}};'.format(n_edges + 1, n_edges + 2))
##geom.add_raw_code('Field[1].FacesList = {rs0:rs7};')
geom.add_raw_code('Field[1].NNodesByEdge = 100;')

# Add new features to refine the mesh around the crystal edges and around melt inclusion
# This will be done by adding raw gmsh code


geom.add_raw_code('Field[2] = Threshold;')
geom.add_raw_code('Field[2].IField = 1;')
geom.add_raw_code('Field[2].DistMax = {0};'.format(2.0*mi_r))
geom.add_raw_code('Field[2].DistMin = {0};'.format(mi_r))
geom.add_raw_code('Field[2].LcMin = {0};'.format(mesh_hres))
geom.add_raw_code('Field[2].LcMax = {0};'.format(mesh_lres))

geom.add_raw_code('Background Field = 2;')

geom.add_raw_code('Mesh.CharacteristicLengthExtendFromBoundary = 0;')
geom.add_raw_code('Mesh.CharacteristicLengthFromPoints = 0;')
geom.add_raw_code('Mesh.CharacteristicLengthFromCurvature = 0;')

code = geom.get_code()

with open(f_mesh + '.geo','w') as out:
    out.write(code)

out.close() 

# Create pandas dataframe of melt inclusion coordinates, radius, tol and crystal size

#pd.DataFrame({'Ol_x': ol_x, 'Ol_y': ol_y, 'Ol_z': ol_z, 'MI_cen_x': mi_x, 'MI_cen_y': mi_y, 'MI_cen_z': mi_z, 'MI_radius': mi_r, 'tolerance': tol}, index = [0]).to_csv(f_coords + '.csv', sep=',')
