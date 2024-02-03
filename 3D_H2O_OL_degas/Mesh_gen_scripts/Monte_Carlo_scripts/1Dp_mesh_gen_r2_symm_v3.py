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

import xtal_geometries as xtal
import olivine_geom as olgeom
import Sphere_geometry as sphg


# script to create comparable 1D meshes for model runs

f_mesh = sys.argv[1]
#f_coords = sys.argv[2]
f_dat = sys.argv[2]
mesh_name = sys.argv[3]
dat_filt = sys.argv[4]

# import dataframe

df_dat = pd.read_csv(f_dat)

# filter dataframe based on mesh name

df_f = df_dat[df_dat['Mesh'] == 'OL_mesh_{0}'.format(mesh_name)]

df_f.to_csv('{0}.csv'.format(dat_filt), sep=',')


# parameters I need

# olivine z
# mi_r
# mi_a dist

# read in parameters from file 

# create coordinates for edge, mi_centre, mi_edge

# add line

# define physical grouping

z = df_f['Ol_z'].values[0] # half length in microns
mi_r = df_f['MI_radius'].values[0] # radius of mi inclusion - microns

dist_a_centre = df_f['MI_centre_a_dist_pve'].values[0]

dist_a_edge = dist_a_centre - mi_r

mi_x = df_f['MI_cen_x'].values[0]
mi_y = df_f['MI_cen_y'].values[0]
mi_z = df_f['MI_cen_z'].values[0]

fc_pve_x = df_f['Face_pve_x'].values[0]
fc_pve_y = df_f['Face_pve_y'].values[0]
fc_pve_z = df_f['Face_pve_z'].values[0]

tol = mi_r/10.0  # edge tolerance
mesh_lres = z/60.0
mesh_hres = mi_r/12.0

pl_th = mesh_lres # plane thickness

geom = pygmsh.built_in.Geometry()

# mi centre
mi_c1 = geom.add_point([mi_x, mi_y + pl_th, mi_z], mesh_hres)
mi_c2 = geom.add_point([mi_x, mi_y - pl_th, mi_z], mesh_hres)

# mi boundary

mi_b1 = geom.add_point([mi_x + mi_r, mi_y + pl_th, mi_z], mesh_hres)
mi_b2 = geom.add_point([mi_x + mi_r, mi_y - pl_th, mi_z], mesh_hres)

# crystal exterior

ol_e1 = geom.add_point([fc_pve_x, fc_pve_y + pl_th, fc_pve_z], mesh_lres)
ol_e2 = geom.add_point([fc_pve_x, fc_pve_y - pl_th, fc_pve_z], mesh_lres)

# Add lines joining points - melt inclusion
mi_l1 = geom.add_line(mi_c1, mi_b1)
mi_l2 = geom.add_line(mi_b1, mi_b2)
mi_l3 = geom.add_line(mi_b2, mi_c2)
mi_l4 = geom.add_line(mi_c2, mi_c1)

# boundary


# Add lines joining points - olivine

ol_l1 = geom.add_line(mi_b2, ol_e2)
ol_l2 = geom.add_line(ol_e2, ol_e1)
ol_l3 = geom.add_line(ol_e1, mi_b1)

mi_ll1 = geom.add_line_loop([mi_l1, mi_l2, mi_l3, mi_l4])
mi_s1 = geom.add_plane_surface(mi_ll1, holes=None)

ol_ll1 = geom.add_line_loop([mi_l2, ol_l1, ol_l2, ol_l3])
ol_s1 = geom.add_plane_surface(ol_ll1, holes=None)


# melt inclusion

mi_ps = geom.add_physical(mi_l2)
mi_pv = geom.add_physical(mi_s1)

ol_ps = geom.add_physical(ol_l2)
ol_pv = geom.add_physical(ol_s1)

# Need to generate coordinate reference frame to position melt inclusion

# define olibine size - scale along z axis (positive direction) - in microns 

# Current olivine dimensions
ol_x = 0.706 # Distance along x (a)
ol_y = 0.46 # Distance along y (b)
ol_z = 0.92 # Distance along z (c)
        
sf = z/ol_z

# Scaled olivine dimensions

ol_x *= sf
ol_y *= sf 
ol_z *= sf 

# Defines the maximum range that the melt inclusion can sit in
max_x = ol_x - (mi_r + tol)
max_y = ol_y - (mi_r + tol)
max_z = ol_z - (mi_r + tol)

min_x = 0.0
min_y = 0.0
min_z = 0.0

     
# Add new features to refine the mesh around the crystal edges and around melt inclusion
# This will be done by adding raw gmsh code

geom.add_raw_code('Field[1] = Attractor;')
geom.add_raw_code('Field[1].EdgesList = {l1};') #mi_l2
#geom.add_raw_code('Field[1].FacesList = {rs0:rs7};')
geom.add_raw_code('Field[1].NNodesByEdge = 100;')
#geom.add_raw_code('Field[1].NodesList = {p1};')


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
