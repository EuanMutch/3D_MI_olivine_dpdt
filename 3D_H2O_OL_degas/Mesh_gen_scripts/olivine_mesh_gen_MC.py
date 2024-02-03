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

f_mesh = sys.argv[1]
f_coords = sys.argv[2]

min_olz = 250.0
max_olz = 750.0

min_mir = 10.0
max_mir = 100.0

z = random.uniform(min_olz, max_olz) # half length in microns
mi_r = random.uniform(min_mir, max_mir) # radius of mi inclusion - microns
tol = mi_r/10.0  # edge tolerance
mesh_lres = z/60.0
mesh_hres = mi_r/12.0

geom = pygmsh.built_in.Geometry()

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

# generate central point o melt inclusion

# Test location of melt inclusion 
# Intersection with crystal faces and in interior of crystal

mi_test = False

# test intersections with faces

ol_gm = olgeom.olivine(z)
ol_faces = ol_gm.faces()

while mi_test == False:

    mi_x = random.uniform(min_x, max_x)
    mi_y = random.uniform(min_y, max_y)
    mi_z = random.uniform(min_z, max_z)

    face_test = np.empty([0])

    for x in ol_faces:
        pl_par = sphg.pl_eq_v1(x[0], x[1], x[2])
        collide_test = sphg.pl_sp_int(pl_par[0], pl_par[1], pl_par[2], pl_par[3], mi_x, mi_y, mi_z, mi_r)
        if collide_test == 0:
            # Test for position of melt inclusion
            c_pos = sphg.pl_pos(pl_par[0], pl_par[1], pl_par[2], pl_par[3], 0.0, 0.0, 0.0)
            sp_pos = sphg.pl_pos(pl_par[0], pl_par[1], pl_par[2], pl_par[3], mi_x, mi_y, mi_z)
            sp_c_sign = np.equal(np.sign(c_pos), np.sign(sp_pos))
            if sp_c_sign == True:
                loc_test = True
            else:
                loc_test = False
        else:
            loc_test = False
        face_test = np.append(face_test, loc_test)  
    
    z_test = np.count_nonzero(face_test)

    if z_test != len(face_test):
       mi_test = False
    else:
       mi_test = True
      

# Run tests

mi_centre = [mi_x, mi_y, mi_z]

mi1 = xtal.sphere(geom, mi_centre, mi_r, mi_r, mi_r, mesh_hres)

mi_s = mi1.surface()

mi_sl = geom.add_surface_loop(mi_s)
mi_ps = geom.add_physical(mi_s)

mi_v = geom.add_volume(mi_sl)
mi_pv = geom.add_physical(mi_v)

ol1 = xtal.olivine(geom, z, mesh_lres)
ol_s = ol1.surfaces()

ol_sl = geom.add_surface_loop(ol_s)
ol_ps = geom.add_physical(ol_s)

ol_v = geom.add_volume(ol_sl, holes= [mi_sl])
ol_pv = geom.add_physical(ol_v)

# Add new features to refine the mesh around the crystal edges and around melt inclusion
# This will be done by adding raw gmsh code

geom.add_raw_code('Field[1] = Attractor;')
geom.add_raw_code('Field[1].EdgesList = {l1:l12};')
geom.add_raw_code('Field[1].FacesList = {rs0:rs7};')
geom.add_raw_code('Field[1].NNodesByEdge = 100;')

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

pd.DataFrame({'Ol_x': ol_x, 'Ol_y': ol_y, 'Ol_z': ol_z, 'MI_cen_x': mi_x, 'MI_cen_y': mi_y, 'MI_cen_z': mi_z, 'MI_radius': mi_r, 'tolerance': tol}, index = [0]).to_csv(f_coords + '.csv', sep=',')
