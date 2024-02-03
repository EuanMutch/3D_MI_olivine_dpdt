import pygmsh
import numpy as np
import meshio
import pandas as pd

#mi_r = 0.1; # Mesh resolution

# Add melt inclusion parameters here


# create class objects for different crystal systems...


# Sphere - melt inclusion

# Olivine

# Feldspars

# Define functions for points, faces, surfaces, volumes etc. 

# Add scaling function in here as well? - yes


class sphere(object):
    
    
    def __init__(self, geom, mi_cp, mi_x, mi_y, mi_z, mi_r):
        self.geom = geom
        self.mi_cp = mi_cp # Melt inclusion central point
        self.mi_x = mi_x # Melt inclusion x axis
        self.mi_y = mi_y # Melt inclusion y axis
        self.mi_z = mi_z # Melt inclusion z axis
        self.mi_r = mi_r # Mesh resolution

    def surface(self):
        mi_p0 = self.geom.add_point(self.mi_cp, self.mi_r)
        mi_p1 = self.geom.add_point([self.mi_cp[0] - self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p2 = self.geom.add_point([self.mi_cp[0] + self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p3 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] - self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p4 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] + self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p5 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] - self.mi_z], self.mi_r)
        mi_p6 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] + self.mi_z], self.mi_r)
        
        #points = [mi_p0, mi_p1, mi_p2, mi_p3, mi_p4, mi_p5, mi_p6]
        
        mi_c1 = self.geom.add_ellipse_arc(mi_p2, mi_p0, mi_p0, mi_p4)
        mi_c2 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p4)
        mi_c3 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p3) 
        mi_c4 = self.geom.add_ellipse_arc(mi_p3, mi_p0, mi_p0, mi_p2) 
        mi_c5 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p2) 
        mi_c6 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p1) 
        mi_c7 = self.geom.add_ellipse_arc(mi_p1, mi_p0, mi_p0, mi_p6) 
        mi_c8 = self.geom.add_ellipse_arc(mi_p6, mi_p0, mi_p0, mi_p3) 
        mi_c9 = self.geom.add_ellipse_arc(mi_p4, mi_p0, mi_p0, mi_p5) 
        mi_c10 = self.geom.add_ellipse_arc(mi_p5, mi_p0, mi_p0, mi_p3) 
        mi_c11 = self.geom.add_ellipse_arc(mi_p3, mi_p0, mi_p0, mi_p6)
        mi_c12 = self.geom.add_ellipse_arc(mi_p6, mi_p0, mi_p0, mi_p4)
        mi_c13 = self.geom.add_ellipse_arc(mi_p2, mi_p0, mi_p0, mi_p6)
        
        #arcs = [mi_c1, mi_c2, mi_c3, mi_c4, mi_c5, mi_c6, mi_c7, mi_c8, mi_c9, mi_c10, mi_c11, mi_c12, mi_c13]
        
        mi_ll1 = self.geom.add_line_loop([mi_c9, mi_c5, mi_c1])
        mi_ll2 = self.geom.add_line_loop([mi_c5, -mi_c4, -mi_c10])
        mi_ll3 = self.geom.add_line_loop([mi_c6, mi_c3, -mi_c10])
        mi_ll4 = self.geom.add_line_loop([mi_c2, mi_c9, mi_c6])
        mi_ll5 = self.geom.add_line_loop([mi_c12, -mi_c2, mi_c7])
        mi_ll6 = self.geom.add_line_loop([mi_c7, mi_c8, -mi_c3])
        mi_ll7 = self.geom.add_line_loop([mi_c12, -mi_c1, mi_c13])
        mi_ll8 = self.geom.add_line_loop([mi_c13, mi_c8, mi_c4])
        
        mi_s1 = self.geom.add_surface(mi_ll1)
        mi_s2 = self.geom.add_surface(mi_ll2)
        mi_s3 = self.geom.add_surface(mi_ll3)
        mi_s4 = self.geom.add_surface(mi_ll4)
        mi_s5 = self.geom.add_surface(mi_ll5)
        mi_s6 = self.geom.add_surface(mi_ll6)
        mi_s7 = self.geom.add_surface(mi_ll7)
        mi_s8 = self.geom.add_surface(mi_ll8)
        
        surfaces = [mi_s1, mi_s2, mi_s3, mi_s4, mi_s5, mi_s6, mi_s7, mi_s8]

        return surfaces
        
    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())



class sphere2(object):
    
    
    def __init__(self, geom, mi_cp, mi_x, mi_y, mi_z, mi_r):
        self.geom = geom
        self.mi_cp = mi_cp # Melt inclusion central point
        self.mi_x = mi_x # Melt inclusion x axis
        self.mi_y = mi_y # Melt inclusion y axis
        self.mi_z = mi_z # Melt inclusion z axis
        self.mi_r = mi_r # Mesh resolution

    def surface(self):
        mi_p0 = self.geom.add_point(self.mi_cp, self.mi_r)
        mi_p1 = self.geom.add_point([self.mi_cp[0] - self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p2 = self.geom.add_point([self.mi_cp[0] + self.mi_x, self.mi_cp[1], self.mi_cp[2]], self.mi_r)
        mi_p3 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] - self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p4 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1] + self.mi_y, self.mi_cp[2]], self.mi_r)
        mi_p5 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] - self.mi_z], self.mi_r)
        mi_p6 = self.geom.add_point([self.mi_cp[0], self.mi_cp[1], self.mi_cp[2] + self.mi_z], self.mi_r)
        
        #points = [mi_p0, mi_p1, mi_p2, mi_p3, mi_p4, mi_p5, mi_p6]
        
        mi_c1 = self.geom.add_circle_arc(mi_p2, mi_p0, mi_p4)
        mi_c2 = self.geom.add_circle_arc(mi_p1, mi_p0, mi_p4)
        mi_c3 = self.geom.add_circle_arc(mi_p1, mi_p0, mi_p3) 
        mi_c4 = self.geom.add_circle_arc(mi_p3, mi_p0, mi_p2) 
        mi_c5 = self.geom.add_circle_arc(mi_p5, mi_p0, mi_p2) 
        mi_c6 = self.geom.add_circle_arc(mi_p5, mi_p0, mi_p1) 
        mi_c7 = self.geom.add_circle_arc(mi_p1, mi_p0, mi_p6) 
        mi_c8 = self.geom.add_circle_arc(mi_p6, mi_p0, mi_p3) 
        mi_c9 = self.geom.add_circle_arc(mi_p4, mi_p0, mi_p5) 
        mi_c10 = self.geom.add_circle_arc(mi_p5, mi_p0, mi_p3) 
        mi_c11 = self.geom.add_circle_arc(mi_p3, mi_p0, mi_p6)
        mi_c12 = self.geom.add_circle_arc(mi_p6, mi_p0, mi_p4)
        mi_c13 = self.geom.add_circle_arc(mi_p2, mi_p0, mi_p6)
        
        #arcs = [mi_c1, mi_c2, mi_c3, mi_c4, mi_c5, mi_c6, mi_c7, mi_c8, mi_c9, mi_c10, mi_c11, mi_c12, mi_c13]
        
        mi_ll1 = self.geom.add_line_loop([mi_c9, mi_c5, mi_c1])
        mi_ll2 = self.geom.add_line_loop([mi_c5, -mi_c4, -mi_c10])
        mi_ll3 = self.geom.add_line_loop([mi_c6, mi_c3, -mi_c10])
        mi_ll4 = self.geom.add_line_loop([mi_c2, mi_c9, mi_c6])
        mi_ll5 = self.geom.add_line_loop([mi_c12, -mi_c2, mi_c7])
        mi_ll6 = self.geom.add_line_loop([mi_c7, mi_c8, -mi_c3])
        mi_ll7 = self.geom.add_line_loop([mi_c12, -mi_c1, mi_c13])
        mi_ll8 = self.geom.add_line_loop([mi_c13, mi_c8, mi_c4])
        
        mi_s1 = self.geom.add_surface(mi_ll1)
        mi_s2 = self.geom.add_surface(mi_ll2)
        mi_s3 = self.geom.add_surface(mi_ll3)
        mi_s4 = self.geom.add_surface(mi_ll4)
        mi_s5 = self.geom.add_surface(mi_ll5)
        mi_s6 = self.geom.add_surface(mi_ll6)
        mi_s7 = self.geom.add_surface(mi_ll7)
        mi_s8 = self.geom.add_surface(mi_ll8)
        
        surfaces = [mi_s1, mi_s2, mi_s3, mi_s4, mi_s5, mi_s6, mi_s7, mi_s8]

        return surfaces
        
    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())


 
class olivine(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Length of z axis of olivine - the rest of the crystal will be scaled accordingly. 
        self.m_r = m_r # Mesh resolution
        
        
        ol_x = 1.0 # Distance along x (a)
        ol_y = 0.46 # Distance along y (b)
        ol_z = 0.92 # Distance along z (c)
        
        self.sf = self.z/ol_z
        
        # define olivine coordinates here
        self.p1 = [0.7061357248994403*self.sf, -0.0000000165845153*self.sf, 0.6865771927157608*self.sf]
        self.p2 = [0.4187895557184330*self.sf, 0.4599999885112464*self.sf, 0.4469191173484135*self.sf]
        self.p3 = [0.6266660391977698*self.sf, 0.1704586163692148*self.sf, 0.7865511922113492*self.sf]
        self.p4 = [0.5645677927000951*self.sf, 0.3036563522524822*self.sf, 0.6303102254605011*self.sf]
        self.p5 = [0.5205871277453623*self.sf, 0.0566915522734537*self.sf, 0.9199999857404487*self.sf]
        self.p6 = [-0.7061357248994403*self.sf, -0.0000000161873903*self.sf, 0.6865771927157608*self.sf]
        self.p7 = [-0.4187895559771333*self.sf, 0.4599999887467703*self.sf, 0.4469191173484135*self.sf]
        self.p8 = [-0.6266660392936343*self.sf, 0.1704586167216467*self.sf, 0.7865511922113492*self.sf]
        self.p9 = [-0.5645677928708690*self.sf, 0.3036563525699906*self.sf, 0.6303102254605011*self.sf]
        self.p10 = [-0.5205871277772453*self.sf, 0.0566915525662277*self.sf, 0.9199999857404486*self.sf]
        self.p11 = [0.4187895559771333*self.sf, -0.4599999887467703*self.sf, 0.4469191173484134*self.sf]
        self.p12 = [0.6266660392936343*self.sf, -0.1704586167216468*self.sf, 0.7865511922113491*self.sf]
        self.p13 = [0.5645677928708690*self.sf, -0.3036563525699906*self.sf, 0.6303102254605010*self.sf]
        self.p14 = [0.5205871277772451*self.sf, -0.0566915525662277*self.sf, 0.9199999857404487*self.sf]
        self.p15 = [-0.4187895557184330*self.sf, -0.4599999885112465*self.sf, 0.4469191173484135*self.sf]
        self.p16 = [-0.6266660391977699*self.sf, -0.1704586163692148*self.sf, 0.7865511922113491*self.sf]
        self.p17 = [-0.5645677927000951*self.sf, -0.3036563522524822*self.sf, 0.6303102254605010*self.sf]
        self.p18 = [-0.5205871277453623*self.sf, -0.0566915522734537*self.sf, 0.9199999857404486*self.sf]
        self.p19 = [0.7061357248994403*self.sf, -0.0000000165845154*self.sf, -0.6865771927157607*self.sf]
        self.p20 = [0.4187895557184330*self.sf, 0.4599999885112465*self.sf, -0.4469191173484135*self.sf]
        self.p21 = [0.6266660391977699*self.sf, 0.1704586163692148*self.sf, -0.7865511922113491*self.sf]
        self.p22 = [0.5645677927000951*self.sf, 0.3036563522524822*self.sf, -0.6303102254605010*self.sf]
        self.p23 = [0.5205871277453623*self.sf, 0.0566915522734537*self.sf, -0.9199999857404486*self.sf]
        self.p24 = [-0.7061357248994404*self.sf, -0.0000000161873904*self.sf, -0.6865771927157607*self.sf]
        self.p25 = [-0.4187895559771333*self.sf, 0.4599999887467703*self.sf, -0.4469191173484134*self.sf]
        self.p26 = [-0.6266660392936343*self.sf, 0.1704586167216468*self.sf, -0.7865511922113491*self.sf]
        self.p27 = [-0.5645677928708690*self.sf, 0.3036563525699906*self.sf, -0.6303102254605010*self.sf]
        self.p28 = [-0.5205871277772451*self.sf, 0.0566915525662277*self.sf, -0.9199999857404487*self.sf]
        self.p29 = [0.4187895559771333*self.sf, -0.4599999887467703*self.sf, -0.4469191173484135*self.sf]
        self.p30 = [0.6266660392936343*self.sf, -0.1704586167216467*self.sf, -0.7865511922113492*self.sf]
        self.p31 = [0.5645677928708690*self.sf, -0.3036563525699906*self.sf, -0.6303102254605011*self.sf]
        self.p32 = [0.5205871277772453*self.sf, -0.0566915525662277*self.sf, -0.9199999857404486*self.sf]
        self.p33 = [-0.4187895557184330*self.sf, -0.4599999885112464*self.sf, -0.4469191173484135*self.sf]
        self.p34 = [-0.6266660391977698*self.sf, -0.1704586163692148*self.sf, -0.7865511922113492*self.sf]
        self.p35 = [-0.5645677927000951*self.sf, -0.3036563522524822*self.sf, -0.6303102254605011*self.sf]
        self.p36 = [-0.5205871277453623*self.sf, -0.0566915522734537*self.sf, -0.9199999857404487*self.sf]

        
        	# Create olivine points
    def surfaces(self):
        ol_p1 = self.geom.add_point(self.p1, self.m_r)
        ol_p2 = self.geom.add_point(self.p2, self.m_r)
        ol_p3 = self.geom.add_point(self.p3, self.m_r)
        ol_p4 = self.geom.add_point(self.p4, self.m_r)
        ol_p5 = self.geom.add_point(self.p5, self.m_r)
        ol_p6 = self.geom.add_point(self.p6, self.m_r)
        ol_p7 = self.geom.add_point(self.p7, self.m_r)
        ol_p8 = self.geom.add_point(self.p8, self.m_r)
        ol_p9 = self.geom.add_point(self.p9, self.m_r)
        ol_p10 = self.geom.add_point(self.p10, self.m_r)
        ol_p11 = self.geom.add_point(self.p11, self.m_r)
        ol_p12 = self.geom.add_point(self.p12, self.m_r)
        ol_p13 = self.geom.add_point(self.p13, self.m_r)
        ol_p14 = self.geom.add_point(self.p14, self.m_r)
        ol_p15 = self.geom.add_point(self.p15, self.m_r)
        ol_p16 = self.geom.add_point(self.p16, self.m_r)
        ol_p17 = self.geom.add_point(self.p17, self.m_r)
        ol_p18 = self.geom.add_point(self.p18, self.m_r)
        ol_p19 = self.geom.add_point(self.p19, self.m_r)
        ol_p20 = self.geom.add_point(self.p20, self.m_r)
        ol_p21 = self.geom.add_point(self.p21, self.m_r)
        ol_p22 = self.geom.add_point(self.p22, self.m_r)
        ol_p23 = self.geom.add_point(self.p23, self.m_r)
        ol_p24 = self.geom.add_point(self.p24, self.m_r)
        ol_p25 = self.geom.add_point(self.p25, self.m_r)
        ol_p26 = self.geom.add_point(self.p26, self.m_r)
        ol_p27 = self.geom.add_point(self.p27, self.m_r)
        ol_p28 = self.geom.add_point(self.p28, self.m_r)
        ol_p29 = self.geom.add_point(self.p29, self.m_r)
        ol_p30 = self.geom.add_point(self.p30, self.m_r)
        ol_p31 = self.geom.add_point(self.p31, self.m_r)
        ol_p32 = self.geom.add_point(self.p32, self.m_r)
        ol_p33 = self.geom.add_point(self.p33, self.m_r)
        ol_p34 = self.geom.add_point(self.p34, self.m_r)
        ol_p35 = self.geom.add_point(self.p35, self.m_r)
        ol_p36 = self.geom.add_point(self.p36, self.m_r)
        
        #return [ol_p1, ol_p2, ol_p3, ol_p4, ol_p5, ol_p6, ol_p7, ol_p8, ol_p9, ol_p10, ol_p11, ol_p12, ol_p13, ol_p14, ol_p15, ol_p16, ol_p17, ol_p18, ol_p19, ol_p20, ol_p21, ol_p22, ol_p23, ol_p24, ol_p25, ol_p26, ol_p27, ol_p28, ol_p29, ol_p30, ol_p31, ol_p32, ol_p33, ol_p34, ol_p35, ol_p36]
        
        ol_l1 = self.geom.add_line(ol_p1, ol_p3)
        ol_l2 = self.geom.add_line(ol_p1, ol_p12)
        ol_l3 = self.geom.add_line(ol_p1, ol_p19)
        ol_l4 = self.geom.add_line(ol_p2, ol_p4)
        ol_l5 = self.geom.add_line(ol_p2, ol_p7)
        ol_l6 = self.geom.add_line(ol_p2, ol_p20)
        ol_l7 = self.geom.add_line(ol_p3, ol_p4)
        ol_l8 = self.geom.add_line(ol_p3, ol_p5)
        ol_l9 = self.geom.add_line(ol_p4, ol_p22)
        ol_l10 = self.geom.add_line(ol_p5, ol_p10)
        ol_l11 = self.geom.add_line(ol_p5, ol_p14)
        ol_l12 = self.geom.add_line(ol_p6, ol_p8)
        ol_l13 = self.geom.add_line(ol_p6, ol_p16)
        ol_l14 = self.geom.add_line(ol_p6, ol_p24)
        ol_l15 = self.geom.add_line(ol_p7, ol_p9)
        ol_l16 = self.geom.add_line(ol_p7, ol_p25)
        ol_l17 = self.geom.add_line(ol_p8, ol_p9)
        ol_l18 = self.geom.add_line(ol_p8, ol_p10)
        ol_l19 = self.geom.add_line(ol_p9, ol_p27)
        ol_l20 = self.geom.add_line(ol_p10, ol_p18)
        ol_l21 = self.geom.add_line(ol_p11, ol_p13)
        ol_l22 = self.geom.add_line(ol_p11, ol_p15)
        ol_l23 = self.geom.add_line(ol_p11, ol_p29)
        ol_l24 = self.geom.add_line(ol_p12, ol_p13)
        ol_l25 = self.geom.add_line(ol_p12, ol_p14)
        ol_l26 = self.geom.add_line(ol_p13, ol_p31)
        ol_l27 = self.geom.add_line(ol_p14, ol_p18)
        ol_l28 = self.geom.add_line(ol_p15, ol_p17)
        ol_l29 = self.geom.add_line(ol_p15, ol_p33)
        ol_l30 = self.geom.add_line(ol_p16, ol_p17)
        ol_l31 = self.geom.add_line(ol_p16, ol_p18)
        ol_l32 = self.geom.add_line(ol_p17, ol_p35)
        ol_l33 = self.geom.add_line(ol_p19, ol_p21)
        ol_l34 = self.geom.add_line(ol_p19, ol_p30)
        ol_l35 = self.geom.add_line(ol_p20, ol_p22)
        ol_l36 = self.geom.add_line(ol_p20, ol_p25)
        ol_l37 = self.geom.add_line(ol_p21, ol_p22)
        ol_l38 = self.geom.add_line(ol_p21, ol_p23)
        ol_l39 = self.geom.add_line(ol_p23, ol_p28)
        ol_l40 = self.geom.add_line(ol_p23, ol_p32)
        ol_l41 = self.geom.add_line(ol_p24, ol_p26)
        ol_l42 = self.geom.add_line(ol_p24, ol_p34)
        ol_l43 = self.geom.add_line(ol_p25, ol_p27)
        ol_l44 = self.geom.add_line(ol_p26, ol_p27)
        ol_l45 = self.geom.add_line(ol_p26, ol_p28)
        ol_l46 = self.geom.add_line(ol_p28, ol_p36)
        ol_l47 = self.geom.add_line(ol_p29, ol_p31)
        ol_l48 = self.geom.add_line(ol_p29, ol_p33)
        ol_l49 = self.geom.add_line(ol_p30, ol_p31)
        ol_l50 = self.geom.add_line(ol_p30, ol_p32)
        ol_l51 = self.geom.add_line(ol_p32, ol_p36)
        ol_l52 = self.geom.add_line(ol_p33, ol_p35)
        ol_l53 = self.geom.add_line(ol_p34, ol_p35)
        ol_l54 = self.geom.add_line(ol_p34, ol_p36)
        
        ol_ll1 = self.geom.add_line_loop([ol_l45, ol_l46, -ol_l54, -ol_l42, ol_l41])
        ol_ll2 = self.geom.add_line_loop([ol_l19, -ol_l44, -ol_l41, -ol_l14, ol_l12, ol_l17])
        ol_ll3 = self.geom.add_line_loop([ol_l42, ol_l53, -ol_l32, -ol_l30, -ol_l13, ol_l14])
        ol_ll4 = self.geom.add_line_loop([ol_l52, -ol_l32, -ol_l28, ol_l29])
        ol_ll5 = self.geom.add_line_loop([ol_l43, -ol_l19, -ol_l15, ol_l16])
        ol_ll6 = self.geom.add_line_loop([ol_l13, ol_l31, -ol_l20, -ol_l18, -ol_l12])
        ol_ll7 = self.geom.add_line_loop([ol_l39, ol_l46, -ol_l51, -ol_l40])
        ol_ll8 = self.geom.add_line_loop([ol_l39, -ol_l45, ol_l44, -ol_l43, -ol_l36, ol_l35, -ol_l37, ol_l38])
        ol_ll9 = self.geom.add_line_loop([ol_l51, -ol_l54, ol_l53, -ol_l52, -ol_l48, ol_l47, -ol_l49, ol_l50])
        ol_ll10 = self.geom.add_line_loop([ol_l40, -ol_l50, -ol_l34, ol_l33, ol_l38])
        ol_ll11 = self.geom.add_line_loop([ol_l49, -ol_l26, -ol_l24, -ol_l2, ol_l3, ol_l34])
        ol_ll12 = self.geom.add_line_loop([ol_l33, ol_l37, -ol_l9, -ol_l7, -ol_l1, ol_l3])
        ol_ll13 = self.geom.add_line_loop([ol_l26, -ol_l47, -ol_l23, ol_l21])
        ol_ll14 = self.geom.add_line_loop([ol_l35, -ol_l9, -ol_l4, ol_l6])
        ol_ll15 = self.geom.add_line_loop([ol_l1, ol_l8, ol_l11, -ol_l25, -ol_l2])
        ol_ll16 = self.geom.add_line_loop([ol_l11, ol_l27, -ol_l20, -ol_l10])
        ol_ll17 = self.geom.add_line_loop([ol_l30, -ol_l28, -ol_l22, ol_l21, -ol_l24, ol_l25, ol_l27, -ol_l31])
        ol_ll18 = self.geom.add_line_loop([ol_l8, ol_l10, -ol_l18, ol_l17, -ol_l15, -ol_l5, ol_l4, -ol_l7])
        ol_ll19 = self.geom.add_line_loop([ol_l29, -ol_l48, -ol_l23, ol_l22])
        ol_ll20 = self.geom.add_line_loop([ol_l16, -ol_l36, -ol_l6, ol_l5])
        
        ol_s1 = self.geom.add_plane_surface(ol_ll1, holes=None)
        ol_s2 = self.geom.add_plane_surface(ol_ll2, holes=None)
        ol_s3 = self.geom.add_plane_surface(ol_ll3, holes=None)
        ol_s4 = self.geom.add_plane_surface(ol_ll4, holes=None)
        ol_s5 = self.geom.add_plane_surface(ol_ll5, holes=None)
        ol_s6 = self.geom.add_plane_surface(ol_ll6, holes=None)
        ol_s7 = self.geom.add_plane_surface(ol_ll7, holes=None)
        ol_s8 = self.geom.add_plane_surface(ol_ll8, holes=None)
        ol_s9 = self.geom.add_plane_surface(ol_ll9, holes=None)
        ol_s10 = self.geom.add_plane_surface(ol_ll10, holes=None)
        ol_s11 = self.geom.add_plane_surface(ol_ll11, holes=None)
        ol_s12 = self.geom.add_plane_surface(ol_ll12, holes=None)
        ol_s13 = self.geom.add_plane_surface(ol_ll13, holes=None)
        ol_s14 = self.geom.add_plane_surface(ol_ll14, holes=None)
        ol_s15 = self.geom.add_plane_surface(ol_ll15, holes=None)
        ol_s16 = self.geom.add_plane_surface(ol_ll16, holes=None)
        ol_s17 = self.geom.add_plane_surface(ol_ll17, holes=None)
        ol_s18 = self.geom.add_plane_surface(ol_ll18, holes=None)
        ol_s19 = self.geom.add_plane_surface(ol_ll19, holes=None)
        ol_s20 = self.geom.add_plane_surface(ol_ll20, holes=None)
        
        
        surfaces = [ol_s1, ol_s2, ol_s3, ol_s4, ol_s5, ol_s6, ol_s7, ol_s8, ol_s9, ol_s10, ol_s11, ol_s12, ol_s13, ol_s14, ol_s15, ol_s16, ol_s17, ol_s18, ol_s19, ol_s20]
        
        return surfaces
        
    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())

    def surface_lines(self, npt):
        list_l1 = [self.p1, self.p1, self.p1, self.p2, self.p2, self.p2, self.p3, self.p3, self.p4, self.p5, self.p5, self.p6, self.p6, self.p6, self.p7, self.p7, self.p8, self.p8, self.p9, self.p10, self.p11, self.p11, self.p11, self.p12, self.p12, self.p13, self.p14, self.p15, self.p15, self.p16, self.p16, self.p17, self.p19, self.p19, self.p20, self.p20, self.p21, self.p21, self.p23, self.p23, self.p24, self.p24, self.p25, self.p26, self.p26, self.p28, self.p29, self.p29, self.p30, self.p30, self.p32, self.p33, self.p34, self.p34]

        list_l2 = [self.p3, self.p12, self.p19, self.p4, self.p7, self.p20, self.p4, self.p5, self.p22, self.p10, self.p14, self.p8, self.p16, self.p24, self.p9, self.p25, self.p9, self.p10, self.p27, self.p18, self.p13, self.p15, self.p29, self.p13, self.p14, self.p31, self.p18, self.p17, self.p33, self.p17, self.p18, self.p35, self.p21, self.p30, self.p22, self.p25, self.p22, self.p23, self.p28, self.p32, self.p26, self.p34, self.p27, self.p27, self.p28, self.p36, self.p31, self.p33, self.p31, self.p32, self.p36, self.p35, self.p35, self.p36]
        
        
        zz = 0    
        while zz < len(list_l1):
            # calculate vector between two points
            vec = [list_l2[zz][0] - list_l1[zz][0], list_l2[zz][1] - list_l1[zz][1], list_l2[zz][2] - list_l1[zz][2]]

            s_pt = 1.0/npt
            s_ptx = s_pt
            # Then calculate coordinates of points along vector
            while s_ptx < 1.0:
                newp = [list_l1[zz][0] + vec[0]*s_ptx, list_l1[zz][1] + vec[1]*s_ptx, list_l1[zz][2] + vec[2]*s_ptx]  
                # Create point geometry
                geom.add_point(newp, self.m_r)

                s_ptx += s_pt
            
            zz += 1   

# plagioclase geometries

# May need to change the coordiantes of fsp_656 to that of orientated file. Need to check. 

class feldspar_656(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Length of z axis of olivine - the rest of the crystal will be scaled accordingly. 
        self.m_r = m_r # Mesh resolution
        
        
        #fsp_x = 1.0 # Distance along x (a)
        fsp_y = 0.45 # Distance along y (b)
        fsp_z = 1.02 # Distance along z (c)
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5543632348405863*self.sf, -0.3377289222132697*self.sf, 0.8007513065339996*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([0.4251258244816719*self.sf, -0.3233829882934972*self.sf, -0.6064248951253192*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.5543632201726041*self.sf, 0.3377289115939922*self.sf, -0.8007513172535427*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([-0.4251258065807629*self.sf, 0.3233829773153498*self.sf, 0.6064249196068604*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([0.6020005884654910*self.sf, 0.1826260440237923*self.sf, 0.9321272745920318*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([0.4272411605936760*self.sf, 0.2020251266666883*self.sf, -0.9707064887421498*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.0936535631128353*self.sf, 0.5590291463518549*self.sf, -0.9712536891890973*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([0.0811058679919058*self.sf, 0.5396300633500893*self.sf, 0.9315801093461680*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([0.4089795335487360*self.sf, 0.3223687943073734*self.sf, 1.0731894449616057*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([0.0993674771359377*self.sf, 0.4192864066875505*self.sf, -1.1123158488391285*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([0.0843502846106453*self.sf, 0.5376403496598434*self.sf, 0.9360172203597096*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([0.4239967260740287*self.sf, 0.2040148513350814*self.sf, -0.9751436242372332*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([0.0936535890721052*self.sf, -0.5590291374991303*self.sf, 0.9712536566574895*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.0811058387997097*self.sf, -0.5396300548562344*self.sf, -0.9315801066766928*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.6020005889086724*self.sf, -0.1826260736910464*self.sf, -0.9321273050914994*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.4272411578039325*self.sf, -0.2020251566928124*self.sf, 0.9707064934437674*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([-0.0993674762859302*self.sf, -0.4192864024489025, 1.1123158278307048*self.sf], self.m_r)*self.sf
        fsp_p18 = self.geom.add_point([-0.4089795382186211*self.sf, -0.3223687981219976*self.sf, -1.0731894655451708*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([-0.4239967410413864*self.sf, -0.2040148701732503*self.sf, 0.9751436044462394*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([-0.0843502734631642*self.sf, -0.5376403303976498*self.sf, -0.9360172421607063*self.sf], self.m_r)
        
        fsp_l21 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l22 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l23 = self.geom.add_line(fsp_p1, fsp_p13)
        fsp_l24 = self.geom.add_line(fsp_p2, fsp_p6)
        fsp_l25 = self.geom.add_line(fsp_p2, fsp_p14)
        fsp_l26 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l27 = self.geom.add_line(fsp_p3, fsp_p7)
        fsp_l28 = self.geom.add_line(fsp_p3, fsp_p15)
        fsp_l29 = self.geom.add_line(fsp_p4, fsp_p8)
        fsp_l30 = self.geom.add_line(fsp_p4, fsp_p16)
        fsp_l31 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l32 = self.geom.add_line(fsp_p5, fsp_p9)
        fsp_l33 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l34 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l35 = self.geom.add_line(fsp_p7, fsp_p10)
        fsp_l36 = self.geom.add_line(fsp_p8, fsp_p11)
        fsp_l37 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l38 = self.geom.add_line(fsp_p9, fsp_p17)
        fsp_l39 = self.geom.add_line(fsp_p10, fsp_p12)
        fsp_l40 = self.geom.add_line(fsp_p10, fsp_p18)
        fsp_l41 = self.geom.add_line(fsp_p11, fsp_p19)
        fsp_l42 = self.geom.add_line(fsp_p12, fsp_p20)
        fsp_l43 = self.geom.add_line(fsp_p13, fsp_p14)
        fsp_l44 = self.geom.add_line(fsp_p13, fsp_p17)
        fsp_l45 = self.geom.add_line(fsp_p14, fsp_p20)
        fsp_l46 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l47 = self.geom.add_line(fsp_p15, fsp_p18)
        fsp_l48 = self.geom.add_line(fsp_p16, fsp_p19)
        fsp_l49 = self.geom.add_line(fsp_p17, fsp_p19)
        fsp_l50 = self.geom.add_line(fsp_p18, fsp_p20)
        
        fsp_ll1 = self.geom.add_line_loop([fsp_l38, fsp_l49, -fsp_l41, -fsp_l37])
        fsp_ll2 = self.geom.add_line_loop([fsp_l44, -fsp_l38, -fsp_l32, -fsp_l22, fsp_l23])
        fsp_ll3 = self.geom.add_line_loop([fsp_l40, fsp_l50, -fsp_l42, -fsp_l39, -fsp_l13])
        fsp_ll4 = self.geom.add_line_loop([fsp_l21, fsp_l25, -fsp_l43, -fsp_l23])
        fsp_ll5 = self.geom.add_line_loop([fsp_l24, -fsp_l31, -fsp_l22, fsp_l21])
        fsp_ll6 = self.geom.add_line_loop([fsp_l27, fsp_l34, -fsp_l29, -fsp_l26])
        fsp_ll7 = self.geom.add_line_loop([fsp_l28, fsp_l46, -fsp_l30, -fsp_l26])
        fsp_ll8 = self.geom.add_line_loop([fsp_l40, -fsp_l47, -fsp_l28, fsp_l27, fsp_l35])
        fsp_ll9 = self.geom.add_line_loop([fsp_l33, -fsp_l39, -fsp_l35, fsp_l34, fsp_l36, -fsp_l37, -fsp_l32, fsp_l31])
        fsp_ll10 = self.geom.add_line_loop([fsp_l41, -fsp_l48, -fsp_l30, fsp_l29, fsp_l36])
        fsp_ll11 = self.geom.add_line_loop([fsp_l49, -fsp_l48, -fsp_l46, fsp_l47, fsp_l50, -fsp_l45, -fsp_l43, fsp_l44])
        fsp_ll12 = self.geom.add_line_loop([fsp_l25, fsp_l45, -fsp_l42, -fsp_l33, -fsp_l24])

        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)

        
        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12]
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())


# Need to change the coordinates of flsp 657 to those from the orientated file.

class feldspar_657(object):
    
    
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Length of z axis of olivine - the rest of the crystal will be scaled accordingly. 
        self.m_r = m_r # Mesh resolution
        
        
        #fsp_x = 1.0 # Distance along x (a)
        fsp_y = 0.45 # Distance along y (b)
        fsp_z = 1.0 # Distance along z (c)
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.0850491638693921*self.sf, 0.5756016323583425*self.sf, 0.8286327720724211*self.sf], self.m_r)
        fsp_p2 = self.geom.add_point([-0.0115193465233267*self.sf, 0.5836470835977820*self.sf, -0.5930940851772852*self.sf], self.m_r)
        fsp_p3 = self.geom.add_point([-0.0850491635551734*self.sf, -0.5756016134159097*self.sf, -0.8286327812328211*self.sf], self.m_r)
        fsp_p4 = self.geom.add_point([0.0115193492331525*self.sf, -0.5836470648549381*self.sf, 0.5930941112861082*self.sf], self.m_r)
        fsp_p5 = self.geom.add_point([-0.3677842383381886*self.sf, 0.3311061239908379*self.sf, 0.9882012761205778*self.sf], self.m_r)
        fsp_p6 = self.geom.add_point([-0.4895352817884513*self.sf, 0.3412496177596306*self.sf, -0.8042746320970990*self.sf], self.m_r)
        fsp_p7 = self.geom.add_point([-0.3747643036728346*self.sf, 0.4032171945846074*self.sf, -0.8447173936000539*self.sf], self.m_r)
        fsp_p8 = self.geom.add_point([-0.5290336376452370*self.sf, -0.2866499494438379*self.sf, -0.9270348203764892*self.sf], self.m_r)
        fsp_p9 = self.geom.add_point([-0.4072825925931848*self.sf, -0.2967934433460785*self.sf, 0.8654411114234576*self.sf], self.m_r)
        fsp_p10 = self.geom.add_point([-0.2947543717454946*self.sf, -0.3700284741979940*self.sf, 0.8903811871216352*self.sf], self.m_r)
        fsp_p11 = self.geom.add_point([-0.3706213509822161*self.sf, 0.1600741416459575*self.sf, 1.0709108970098904*self.sf], self.m_r)
        fsp_p12 = self.geom.add_point([-0.5261965254214045*self.sf, -0.1156179924301354*self.sf, -1.0097444290158593*self.sf], self.m_r)
        fsp_p13 = self.geom.add_point([-0.2006920515361074*self.sf, 0.0494816396599514*self.sf, 1.1085730013211323*self.sf], self.m_r)
        fsp_p14 = self.geom.add_point([-0.3528804026530329*self.sf, -0.0220405015958488*self.sf, -1.0708172083266450*self.sf], self.m_r)
        fsp_p15 = self.geom.add_point([-0.2318481768115107*self.sf, -0.3387819075305624*self.sf, 0.9339569190844070*self.sf], self.m_r)
        fsp_p16 = self.geom.add_point([-0.3217242773776302*self.sf, 0.3662230455946638*self.sf, -0.8962011260899205*self.sf], self.m_r)
        fsp_p17 = self.geom.add_point([0.5290336132547900*self.sf, 0.2866499696586574*self.sf, 0.9270347919682673*self.sf], self.m_r)*self.sf
        fsp_p18 = self.geom.add_point([0.4072825698045237*self.sf, 0.2967934634274501*self.sf, -0.8654411162494090*self.sf], self.m_r)
        fsp_p19 = self.geom.add_point([0.2947543547683388*self.sf, 0.3700284876971762*self.sf, -0.8903811840101974*self.sf], self.m_r)
        fsp_p20 = self.geom.add_point([0.3677842605265381*self.sf, -0.3311061060849770*self.sf, -0.9882013077056566*self.sf], self.m_r)
        fsp_p21 = self.geom.add_point([0.4895353055785912*self.sf, -0.3412495999872183*self.sf, 0.8042746240942899*self.sf], self.m_r)
        fsp_p22 = self.geom.add_point([0.3747643221606282*self.sf, -0.4032171828434580*self.sf, 0.8447173942926415*self.sf], self.m_r)
        fsp_p23 = self.geom.add_point([0.5261965132495373*self.sf, 0.1156179866872595*self.sf, 1.0097444119955694*self.sf], self.m_r)
        fsp_p24 = self.geom.add_point([0.3706213601115924*self.sf, -0.1600741484447564*self.sf, -1.0709109154830110*self.sf], self.m_r)
        fsp_p25 = self.geom.add_point([0.3528804094132100*self.sf, 0.0220405012902180*self.sf, 1.0708171949445484*self.sf], self.m_r)
        fsp_p26 = self.geom.add_point([0.2006920430287449*self.sf, -0.0494816392087548*self.sf, -1.1085730136619267*self.sf], self.m_r)*self.sf
        fsp_p27 = self.geom.add_point([0.3217243129813364*self.sf, -0.3662230473300983*self.sf, 0.8962011107405788*self.sf], self.m_r)
        fsp_p28 = self.geom.add_point([0.2318481394606212*self.sf, 0.3387819094115611*self.sf, -0.9339569294579589*self.sf], self.m_r)

        
        fsp_l29 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l30 = self.geom.add_line(fsp_p1, fsp_p5)
        fsp_l31 = self.geom.add_line(fsp_p1, fsp_p17)
        fsp_l32 = self.geom.add_line(fsp_p2, fsp_p7)
        fsp_l33 = self.geom.add_line(fsp_p2, fsp_p19)
        fsp_l34 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l35 = self.geom.add_line(fsp_p3, fsp_p8)
        fsp_l36 = self.geom.add_line(fsp_p3, fsp_p20)
        fsp_l37 = self.geom.add_line(fsp_p4, fsp_p10)
        fsp_l38 = self.geom.add_line(fsp_p4, fsp_p22)
        fsp_l39 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l40 = self.geom.add_line(fsp_p5, fsp_p11)
        fsp_l41 = self.geom.add_line(fsp_p6, fsp_p7)
        fsp_l42 = self.geom.add_line(fsp_p6, fsp_p12)
        fsp_l43 = self.geom.add_line(fsp_p7, fsp_p16)
        fsp_l44 = self.geom.add_line(fsp_p8, fsp_p9)
        fsp_l45 = self.geom.add_line(fsp_p8, fsp_p12)
        fsp_l46 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l47 = self.geom.add_line(fsp_p9, fsp_p11)
        fsp_l48 = self.geom.add_line(fsp_p10, fsp_p15)
        fsp_l49 = self.geom.add_line(fsp_p11, fsp_p13)
        fsp_l50 = self.geom.add_line(fsp_p12, fsp_p14)
        fsp_l51 = self.geom.add_line(fsp_p13, fsp_p15)
        fsp_l52 = self.geom.add_line(fsp_p13, fsp_p25)
        fsp_l53 = self.geom.add_line(fsp_p14, fsp_p16)
        fsp_l54 = self.geom.add_line(fsp_p14, fsp_p26)
        fsp_l55 = self.geom.add_line(fsp_p15, fsp_p27)
        fsp_l56 = self.geom.add_line(fsp_p16, fsp_p28)
        fsp_l57 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l58 = self.geom.add_line(fsp_p17, fsp_p23)
        fsp_l59 = self.geom.add_line(fsp_p18, fsp_p19)
        fsp_l60 = self.geom.add_line(fsp_p18, fsp_p24)
        fsp_l61 = self.geom.add_line(fsp_p19, fsp_p28)
        fsp_l62 = self.geom.add_line(fsp_p20, fsp_p21)
        fsp_l63 = self.geom.add_line(fsp_p20, fsp_p24)
        fsp_l64 = self.geom.add_line(fsp_p21, fsp_p22)
        fsp_l65 = self.geom.add_line(fsp_p21, fsp_p23)
        fsp_l66 = self.geom.add_line(fsp_p22, fsp_p27)
        fsp_l67 = self.geom.add_line(fsp_p23, fsp_p25)
        fsp_l68 = self.geom.add_line(fsp_p24, fsp_p26)
        fsp_l69 = self.geom.add_line(fsp_p25, fsp_p27)
        fsp_l70 = self.geom.add_line(fsp_p26, fsp_p28)

        
        fsp_ll1 = self.geom.add_line_loop([fsp_l55, -fsp_l66, -fsp_l38, fsp_l37, fsp_l48])
        fsp_ll2 = self.geom.add_line_loop([fsp_l51, fsp_l55, -fsp_l69, -fsp_l52])
        fsp_ll3 = self.geom.add_line_loop([fsp_l69, -fsp_l66, -fsp_l64, fsp_l65, fsp_l67])
        fsp_ll4 = self.geom.add_line_loop([fsp_l52, -fsp_l67, -fsp_l58, -fsp_l31, fsp_l30, fsp_l40, fsp_l49])
        fsp_ll5 = self.geom.add_line_loop([fsp_l51, -fsp_l48, -fsp_l46, fsp_l47, fsp_l49])
        fsp_ll6 = self.geom.add_line_loop([fsp_l37, -fsp_l46, -fsp_l44, -fsp_l35, fsp_l34])
        fsp_ll7 = self.geom.add_line_loop([fsp_l47, -fsp_l40, fsp_l39, fsp_l42, -fsp_l45, fsp_l44])
        fsp_ll8 = self.geom.add_line_loop([fsp_l30, fsp_l39, fsp_l41, -fsp_l32, -fsp_l29])
        fsp_ll9 = self.geom.add_line_loop([fsp_l31, fsp_l57, fsp_l59, -fsp_l33, -fsp_l29])
        fsp_ll10 = self.geom.add_line_loop([fsp_l58, -fsp_l65, -fsp_l62, fsp_l63, -fsp_l60, -fsp_l57])
        fsp_ll11 = self.geom.add_line_loop([fsp_l59, fsp_l61, -fsp_l70, -fsp_l68, -fsp_l60])
        fsp_ll12 = self.geom.add_line_loop([fsp_l33, fsp_l61, -fsp_l56, -fsp_l43, -fsp_l32])
        fsp_ll13 = self.geom.add_line_loop([fsp_l56, -fsp_l70, -fsp_l54, fsp_l53])
        fsp_ll14 = self.geom.add_line_loop([fsp_l42, fsp_l50, fsp_l53, -fsp_l43, -fsp_l41])
        fsp_ll15 = self.geom.add_line_loop([fsp_l54, -fsp_l68, -fsp_l63, -fsp_l36, fsp_l35, fsp_l45, fsp_l50])
        fsp_ll16 = self.geom.add_line_loop([fsp_l38, -fsp_l64, -fsp_l62, -fsp_l36, fsp_l34])

        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)
        fsp_s13 = self.geom.add_plane_surface(fsp_ll13, holes=None)
        fsp_s14 = self.geom.add_plane_surface(fsp_ll14, holes=None)
        fsp_s15 = self.geom.add_plane_surface(fsp_ll15, holes=None)
        fsp_s16 = self.geom.add_plane_surface(fsp_ll16, holes=None)

        
        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16] 
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())

class anorthite_743(object):
    
# Need to change the coordiantes to that from orientated file.     
    def __init__(self, geom, z, m_r):
        self.geom = geom
        self.z = z # Length of z axis of olivine - the rest of the crystal will be scaled accordingly. 
        self.m_r = m_r # Mesh resolution
        
        
        fsp_x = 0.56 # Distance along x (a)
        fsp_y = 0.5 # Distance along y (b)
        fsp_z = 0.48 # Distance along z (c)
        
        self.sf = self.z/fsp_z


    def surfaces(self):
        fsp_p1 = self.geom.add_point([0.5342785389845653, -0.1546832235199951, 0.1224330536660508], self.m_r)
        fsp_p2 = self.geom.add_point([0.5391106330277696, -0.1879434835582076, -0.0776893812182014], self.m_r)
        fsp_p3 = self.geom.add_point([0.4943787992470607, -0.2411445167791730, 0.1412478934939521], self.m_r)
        fsp_p4 = self.geom.add_point([0.4992108932902646, -0.2744047768173857, -0.0588745413902999], self.m_r)
        fsp_p5 = self.geom.add_point([-0.5342785389845653, 0.1546832235199951, -0.1224330536660508], self.m_r)
        fsp_p6 = self.geom.add_point([-0.5391106330277696, 0.1879434835582076, 0.0776893812182014], self.m_r)
        fsp_p7 = self.geom.add_point([-0.4943787992470607, 0.2411445167791730, -0.1412478934939521], self.m_r)
        fsp_p8 = self.geom.add_point([-0.4992108932902646, 0.2744047768173857, 0.0588745413902999], self.m_r)
        fsp_p9 = self.geom.add_point([0.5300756432361228, -0.1129712626545864, 0.1814093674606908], self.m_r)
        fsp_p10 = self.geom.add_point([0.4991611590483267, 0.2411938547722703, 0.1888790424749543], self.m_r)
        fsp_p11 = self.geom.add_point([0.5089514630407531, 0.1738052539441853, -0.2165889504290384], self.m_r)
        fsp_p12 = self.geom.add_point([0.5294246959184634, -0.0320374430888905, -0.4799756269307332], self.m_r)
        fsp_p13 = self.geom.add_point([0.5276916440312652, -0.0121831080979484, -0.4795568803554521], self.m_r)
        fsp_p14 = self.geom.add_point([-0.5300756432361228, 0.1129712626545864, -0.1814093674606908], self.m_r)
        fsp_p15 = self.geom.add_point([-0.4991611590483267, -0.2411938547722703, -0.1888790424749543], self.m_r)
        fsp_p16 = self.geom.add_point([-0.5089514630407531, -0.1738052539441853, 0.2165889504290384], self.m_r)
        fsp_p17 = self.geom.add_point([-0.5294246959184634, 0.0320374430888905, 0.4799756269307332], self.m_r)*self.sf
        fsp_p18 = self.geom.add_point([-0.5276916440312652, 0.0121831080979484, 0.4795568803554521], self.m_r)
        fsp_p19 = self.geom.add_point([-0.4452425090381870, 0.2527010664851914, -0.2708717518487827], self.m_r)
        fsp_p20 = self.geom.add_point([-0.0875567124977610, 0.4807753669359448, -0.4303071369751568], self.m_r)
        fsp_p21 = self.geom.add_point([-0.0998266387102427, 0.5652316984617675, 0.0778550567936375], self.m_r)
        fsp_p22 = self.geom.add_point([-0.3850271562689780, 0.3986727222558797, 0.2883047135876014], self.m_r)
        fsp_p23 = self.geom.add_point([-0.1378340930018188, 0.5562925817141924, 0.1781205194603585], self.m_r)
        fsp_p24 = self.geom.add_point([0.4452425090381870, -0.2527010664851914, 0.2708717518487827], self.m_r)
        fsp_p25 = self.geom.add_point([0.0875567124977610, -0.4807753669359448, 0.4303071369751568], self.m_r)
        fsp_p26 = self.geom.add_point([0.0998266387102427, -0.5652316984617675, -0.0778550567936375], self.m_r)*self.sf
        fsp_p27 = self.geom.add_point([0.3850271562689780, -0.3986727222558797, -0.2883047135876014], self.m_r)
        fsp_p28 = self.geom.add_point([0.1378340930018188, -0.5562925817141924, -0.1781205194603585], self.m_r)
        fsp_p29 = self.geom.add_point([0.4036992675088759, 0.3798821627699073, 0.1341251802267428], self.m_r)
        fsp_p30 = self.geom.add_point([0.4157168084107762, 0.2971630493082779, -0.3635844100888118], self.m_r)
        fsp_p31 = self.geom.add_point([0.0248629077209667, 0.4981576005700549, -0.3722117235781828], self.m_r)
        fsp_p32 = self.geom.add_point([0.0121247871714059, 0.5858366064189519, 0.1553408607537549], self.m_r)
        fsp_p33 = self.geom.add_point([0.3780341180594647, 0.2922652290970223, -0.5223589948293154], self.m_r)
        fsp_p34 = self.geom.add_point([0.0434399026859767, 0.5899067989730189, 0.2872859107375331], self.m_r)
        fsp_p35 = self.geom.add_point([-0.4036992675088759, -0.3798821627699073, -0.1341251802267428], self.m_r)
        fsp_p36 = self.geom.add_point([-0.4157168084107762, -0.2971630493082779, 0.3635844100888118], self.m_r)
        fsp_p37 = self.geom.add_point([-0.0248629077209667, -0.4981576005700549, 0.3722117235781828], self.m_r)
        fsp_p38 = self.geom.add_point([-0.0121247871714059, -0.5858366064189519, -0.1553408607537549], self.m_r)
        fsp_p39 = self.geom.add_point([-0.3780341180594647, -0.2922652290970223, 0.5223589948293154], self.m_r)
        fsp_p40 = self.geom.add_point([-0.0434399026859767, -0.5899067989730189, -0.2872859107375331], self.m_r)
        fsp_p41 = self.geom.add_point([0.4784946473916337, -0.1251027919138948, 0.3174824749485568], self.m_r)
        fsp_p42 = self.geom.add_point([0.4406490419496569, -0.2071128504309959, 0.3353286816974050], self.m_r)
        fsp_p43 = self.geom.add_point([0.4402892810142589, 0.3125887441047152, 0.3267138005921055], self.m_r)
        fsp_p44 = self.geom.add_point([0.1198065891922738, -0.4116944096522709, 0.4783414613912789], self.m_r)
        fsp_p45 = self.geom.add_point([0.0688122729214938, 0.5291530616509266, 0.4846436019601287], self.m_r)*self.sf
        fsp_p46 = self.geom.add_point([-0.3953985360093574, -0.1113392223731165, 0.6973758600600728], self.m_r)
        fsp_p47 = self.geom.add_point([-0.2780677435445308, 0.3079689445821555, 0.6392624030218961], self.m_r)
        fsp_p48 = self.geom.add_point([-0.4075143384628512, 0.0274628541097272, 0.7003033262378963], self.m_r)
        fsp_p49 = self.geom.add_point([-0.4784946473916337, 0.1251027919138948, -0.3174824749485568], self.m_r)
        fsp_p50 = self.geom.add_point([-0.4406490419496569, 0.2071128504309959, -0.3353286816974050], self.m_r)
        fsp_p51 = self.geom.add_point([-0.4402892810142589, -0.3125887441047152, -0.3267138005921055], self.m_r)
        fsp_p52 = self.geom.add_point([-0.1198065891922738, 0.4116944096522709, -0.4783414613912789], self.m_r)
        fsp_p53 = self.geom.add_point([-0.0688122729214938, -0.5291530616509266, -0.4846436019601287], self.m_r)
        fsp_p54 = self.geom.add_point([0.3953985360093574, 0.1113392223731165, -0.6973758600600728], self.m_r)*self.sf
        fsp_p55 = self.geom.add_point([0.2780677435445308, -0.3079689445821555, -0.6392624030218961], self.m_r)
        fsp_p56 = self.geom.add_point([0.4075143384628512, -0.0274628541097272, -0.7003033262378963], self.m_r)
        fsp_p57 = self.geom.add_point([0.4575607137126042, 0.3166508776036849, 0.2811509439710933], self.m_r)
        fsp_p58 = self.geom.add_point([0.4171768940053555, 0.3700822866354242, 0.2174810761065769], self.m_r)
        fsp_p59 = self.geom.add_point([-0.4575607137126042, -0.3166508776036849, -0.2811509439710933], self.m_r)
        fsp_p60 = self.geom.add_point([-0.4171768940053555, -0.3700822866354242, -0.2174810761065769], self.m_r)
        fsp_p61 = self.geom.add_point([-0.0785764300032806, 0.4791193612964130, -0.4438059903018485], self.m_r)
        fsp_p62 = self.geom.add_point([0.0785764300032806, -0.4791193612964130, 0.4438059903018485], self.m_r)
        fsp_p63 = self.geom.add_point([0.0762052460250527, -0.4555863642333433, 0.4770791664011019], self.m_r)
        fsp_p64 = self.geom.add_point([-0.0762052460250527, 0.4555863642333433, -0.4770791664011019], self.m_r)
        fsp_p65 = self.geom.add_point([0.0659747147282316, 0.5748267489276897, 0.3667912331416582], self.m_r)
        fsp_p66 = self.geom.add_point([-0.0659747147282316, -0.5748267489276897, -0.3667912331416582], self.m_r)
        fsp_p67 = self.geom.add_point([-0.3888986575002198, -0.1844392750173213, 0.6748135227650793], self.m_r)
        fsp_p68 = self.geom.add_point([0.3888986575002198, 0.1844392750173213, -0.6748135227650793], self.m_r)
        fsp_p69 = self.geom.add_point([0.0433070517712530, 0.5896320107202505, 0.3034951372685833], self.m_r)
        fsp_p70 = self.geom.add_point([-0.0433070517712530, -0.5896320107202505, -0.3034951372685833], self.m_r)
        fsp_p71 = self.geom.add_point([-0.3978709102724138, -0.1595813455013852, 0.6842346207084005], self.m_r)
        fsp_p72 = self.geom.add_point([0.3978709102724138, 0.1595813455013852, -0.6842346207084005], self.m_r)
        fsp_p73 = self.geom.add_point([-0.3460584059033321, 0.2954551125445693, 0.5922035014209278], self.m_r)*self.sf
        fsp_p74 = self.geom.add_point([-0.3791190166711064, 0.3035745773010317, 0.5336876095482047], self.m_r)
        fsp_p75 = self.geom.add_point([0.3460584059033321, -0.2954551125445693, -0.5922035014209278], self.m_r)
        fsp_p76 = self.geom.add_point([0.3791190166711064, -0.3035745773010317, -0.5336876095482047], self.m_r)
        fsp_p77 = self.geom.add_point([-0.4461439697261558, 0.0785733147778230, 0.6393991439188552], self.m_r)
        fsp_p78 = self.geom.add_point([-0.4792045804939291, 0.0866927795342853, 0.5808832520461336], self.m_r)
        fsp_p79 = self.geom.add_point([0.4461439697261558, -0.0785733147778230, -0.6393991439188552], self.m_r)
        fsp_p80 = self.geom.add_point([0.4792045804939291, -0.0866927795342853, -0.5808832520461336], self.m_r)


        
        fsp_l81 = self.geom.add_line(fsp_p1, fsp_p2)
        fsp_l82 = self.geom.add_line(fsp_p1, fsp_p3)
        fsp_l83 = self.geom.add_line(fsp_p1, fsp_p9)
        fsp_l84 = self.geom.add_line(fsp_p2, fsp_p4)
        fsp_l85 = self.geom.add_line(fsp_p2, fsp_p12)
        fsp_l86 = self.geom.add_line(fsp_p3, fsp_p4)
        fsp_l87 = self.geom.add_line(fsp_p3, fsp_p24)
        fsp_l88 = self.geom.add_line(fsp_p4, fsp_p27)
        fsp_l89 = self.geom.add_line(fsp_p5, fsp_p6)
        fsp_l90 = self.geom.add_line(fsp_p5, fsp_p7)
        fsp_l91 = self.geom.add_line(fsp_p5, fsp_p14)
        fsp_l92 = self.geom.add_line(fsp_p6, fsp_p8)
        fsp_l93 = self.geom.add_line(fsp_p6, fsp_p17)
        fsp_l94 = self.geom.add_line(fsp_p7, fsp_p8)
        fsp_l95 = self.geom.add_line(fsp_p7, fsp_p19)
        fsp_l96 = self.geom.add_line(fsp_p8, fsp_p22)
        fsp_l97 = self.geom.add_line(fsp_p9, fsp_p10)
        fsp_l98 = self.geom.add_line(fsp_p9, fsp_p41)
        fsp_l99 = self.geom.add_line(fsp_p10, fsp_p11)
        fsp_l100 = self.geom.add_line(fsp_p10, fsp_p57)
        fsp_l101 = self.geom.add_line(fsp_p11, fsp_p13)
        fsp_l102 = self.geom.add_line(fsp_p11, fsp_p30)
        fsp_l103 = self.geom.add_line(fsp_p12, fsp_p13)
        fsp_l104 = self.geom.add_line(fsp_p12, fsp_p80)
        fsp_l105 = self.geom.add_line(fsp_p13, fsp_p72)
        fsp_l106 = self.geom.add_line(fsp_p14, fsp_p15)
        fsp_l107 = self.geom.add_line(fsp_p14, fsp_p49)
        fsp_l108 = self.geom.add_line(fsp_p15, fsp_p16)
        fsp_l109 = self.geom.add_line(fsp_p15, fsp_p59)
        fsp_l110 = self.geom.add_line(fsp_p16, fsp_p18)
        fsp_l111 = self.geom.add_line(fsp_p16, fsp_p36)
        fsp_l112 = self.geom.add_line(fsp_p17, fsp_p18)
        fsp_l113 = self.geom.add_line(fsp_p17, fsp_p78)
        fsp_l114 = self.geom.add_line(fsp_p18, fsp_p71)
        fsp_l115 = self.geom.add_line(fsp_p19, fsp_p20)
        fsp_l116 = self.geom.add_line(fsp_p19, fsp_p50)
        fsp_l117 = self.geom.add_line(fsp_p20, fsp_p21)
        fsp_l118 = self.geom.add_line(fsp_p20, fsp_p61)
        fsp_l119 = self.geom.add_line(fsp_p21, fsp_p23)
        fsp_l120 = self.geom.add_line(fsp_p21, fsp_p32)
        fsp_l121 = self.geom.add_line(fsp_p22, fsp_p23)
        fsp_l122 = self.geom.add_line(fsp_p22, fsp_p74)
        fsp_l123 = self.geom.add_line(fsp_p23, fsp_p69)
        fsp_l124 = self.geom.add_line(fsp_p24, fsp_p25)
        fsp_l125 = self.geom.add_line(fsp_p24, fsp_p42)
        fsp_l126 = self.geom.add_line(fsp_p25, fsp_p26)
        fsp_l127 = self.geom.add_line(fsp_p25, fsp_p62)
        fsp_l128 = self.geom.add_line(fsp_p26, fsp_p28)
        fsp_l129 = self.geom.add_line(fsp_p26, fsp_p38)
        fsp_l130 = self.geom.add_line(fsp_p27, fsp_p28)
        fsp_l131 = self.geom.add_line(fsp_p27, fsp_p76)
        fsp_l132 = self.geom.add_line(fsp_p28, fsp_p70)
        fsp_l133 = self.geom.add_line(fsp_p29, fsp_p30)
        fsp_l134 = self.geom.add_line(fsp_p29, fsp_p34)
        fsp_l135 = self.geom.add_line(fsp_p29, fsp_p58)
        fsp_l136 = self.geom.add_line(fsp_p30, fsp_p33)
        fsp_l137 = self.geom.add_line(fsp_p31, fsp_p32)
        fsp_l138 = self.geom.add_line(fsp_p31, fsp_p33)
        fsp_l139 = self.geom.add_line(fsp_p31, fsp_p61)
        fsp_l140 = self.geom.add_line(fsp_p32, fsp_p34)
        fsp_l141 = self.geom.add_line(fsp_p33, fsp_p68)
        fsp_l142 = self.geom.add_line(fsp_p34, fsp_p69)
        fsp_l143 = self.geom.add_line(fsp_p35, fsp_p36)
        fsp_l144 = self.geom.add_line(fsp_p35, fsp_p40)
        fsp_l145 = self.geom.add_line(fsp_p35, fsp_p60)
        fsp_l146 = self.geom.add_line(fsp_p36, fsp_p39)
        fsp_l147 = self.geom.add_line(fsp_p37, fsp_p38)
        fsp_l148 = self.geom.add_line(fsp_p37, fsp_p39)
        fsp_l149 = self.geom.add_line(fsp_p37, fsp_p62)
        fsp_l150 = self.geom.add_line(fsp_p38, fsp_p40)
        fsp_l151 = self.geom.add_line(fsp_p39, fsp_p67)
        fsp_l152 = self.geom.add_line(fsp_p40, fsp_p70)
        fsp_l153 = self.geom.add_line(fsp_p41, fsp_p42)
        fsp_l154 = self.geom.add_line(fsp_p41, fsp_p43)
        fsp_l155 = self.geom.add_line(fsp_p42, fsp_p44)
        fsp_l156 = self.geom.add_line(fsp_p43, fsp_p45)
        fsp_l157 = self.geom.add_line(fsp_p43, fsp_p57)
        fsp_l158 = self.geom.add_line(fsp_p44, fsp_p46)
        fsp_l159 = self.geom.add_line(fsp_p44, fsp_p63)
        fsp_l160 = self.geom.add_line(fsp_p45, fsp_p47)
        fsp_l161 = self.geom.add_line(fsp_p45, fsp_p65)
        fsp_l162 = self.geom.add_line(fsp_p46, fsp_p48)
        fsp_l163 = self.geom.add_line(fsp_p46, fsp_p71)
        fsp_l164 = self.geom.add_line(fsp_p47, fsp_p48)
        fsp_l165 = self.geom.add_line(fsp_p47, fsp_p73)
        fsp_l166 = self.geom.add_line(fsp_p48, fsp_p77)
        fsp_l167 = self.geom.add_line(fsp_p49, fsp_p50)
        fsp_l168 = self.geom.add_line(fsp_p49, fsp_p51)
        fsp_l169 = self.geom.add_line(fsp_p50, fsp_p52)
        fsp_l170 = self.geom.add_line(fsp_p51, fsp_p53)
        fsp_l171 = self.geom.add_line(fsp_p51, fsp_p59)
        fsp_l172 = self.geom.add_line(fsp_p52, fsp_p54)
        fsp_l173 = self.geom.add_line(fsp_p52, fsp_p64)
        fsp_l174 = self.geom.add_line(fsp_p53, fsp_p55)
        fsp_l175 = self.geom.add_line(fsp_p53, fsp_p66)
        fsp_l176 = self.geom.add_line(fsp_p54, fsp_p56)
        fsp_l177 = self.geom.add_line(fsp_p54, fsp_p72)
        fsp_l178 = self.geom.add_line(fsp_p55, fsp_p56)
        fsp_l179 = self.geom.add_line(fsp_p55, fsp_p75)
        fsp_l180 = self.geom.add_line(fsp_p56, fsp_p79)
        fsp_l181 = self.geom.add_line(fsp_p57, fsp_p58)
        fsp_l182 = self.geom.add_line(fsp_p57, fsp_p58)
        fsp_l183 = self.geom.add_line(fsp_p59, fsp_p60)
        fsp_l184 = self.geom.add_line(fsp_p60, fsp_p66)
        fsp_l185 = self.geom.add_line(fsp_p61, fsp_p64)
        fsp_l186 = self.geom.add_line(fsp_p62, fsp_p63)
        fsp_l187 = self.geom.add_line(fsp_p63, fsp_p67)
        fsp_l188 = self.geom.add_line(fsp_p64, fsp_p68)
        fsp_l189 = self.geom.add_line(fsp_p65, fsp_p69)
        fsp_l190 = self.geom.add_line(fsp_p66, fsp_p70)
        fsp_l191 = self.geom.add_line(fsp_p67, fsp_p71)
        fsp_l192 = self.geom.add_line(fsp_p68, fsp_p72)
        fsp_l193 = self.geom.add_line(fsp_p73, fsp_p74)
        fsp_l194 = self.geom.add_line(fsp_p73, fsp_p77)
        fsp_l195 = self.geom.add_line(fsp_p74, fsp_p78)
        fsp_l196 = self.geom.add_line(fsp_p75, fsp_p76)
        fsp_l197 = self.geom.add_line(fsp_p75, fsp_p79)
        fsp_l198 = self.geom.add_line(fsp_p76, fsp_p80)
        fsp_l199 = self.geom.add_line(fsp_p77, fsp_p78)
        fsp_l200 = self.geom.add_line(fsp_p79, fsp_p80)
        
        fsp_ll1 = self.geom.add_line_loop([fsp_l168, -fsp_l171, -fsp_l109, -fsp_l106, fsp_l107])
        fsp_ll2 = self.geom.add_line_loop([fsp_l107, fsp_l167, -fsp_l116, -fsp_l95, -fsp_l90, fsp_l91])
        fsp_ll3 = self.geom.add_line_loop([fsp_l89, fsp_l92, -fsp_l94, -fsp_l90])
        fsp_ll4 = self.geom.add_line_loop([fsp_l93, fsp_l113, -fsp_l195, -fsp_l122, -fsp_l96, -fsp_l92])
        fsp_ll5 = self.geom.add_line_loop([fsp_l93, fsp_l112, -fsp_l110, -fsp_l108, -fsp_l106, -fsp_l91, fsp_l89])
        fsp_ll6 = self.geom.add_line_loop([fsp_l170, fsp_l174, fsp_l178, -fsp_l176, -fsp_l172, -fsp_l169, -fsp_l167, fsp_l168])
        fsp_ll7 = self.geom.add_line_loop([fsp_l170, fsp_l175, -fsp_l184, -fsp_l183, -fsp_l171])
        fsp_ll8 = self.geom.add_line_loop([fsp_l145, fsp_l184, fsp_l190, -fsp_l152, -fsp_l144])
        fsp_ll9 = self.geom.add_line_loop([fsp_l143, -fsp_l111, -fsp_l108, fsp_l109, fsp_l183, -fsp_l145])
        fsp_ll10 = self.geom.add_line_loop([fsp_l110, fsp_l114, -fsp_l191, -fsp_l151, -fsp_l146, -fsp_l111])
        fsp_ll11 = self.geom.add_line_loop([fsp_l95, fsp_l115, fsp_l117, fsp_l119, -fsp_l121, -fsp_l96, -fsp_l94])
        fsp_ll12 = self.geom.add_line_loop([fsp_l115, fsp_l118, fsp_l185, -fsp_l173, -fsp_l169, -fsp_l116])
        fsp_ll13 = self.geom.add_line_loop([fsp_l172, fsp_l177, -fsp_l192, -fsp_l188, -fsp_l173])
        fsp_ll14 = self.geom.add_line_loop([fsp_l185, fsp_l188, -fsp_l141, -fsp_l138, fsp_l139])
        fsp_ll15 = self.geom.add_line_loop([fsp_l141, fsp_l192, -fsp_l105, -fsp_l101, fsp_l102, fsp_l136])
        fsp_ll16 = self.geom.add_line_loop([fsp_l138, -fsp_l136, -fsp_l133, fsp_l134, -fsp_l140, -fsp_l137])
        fsp_ll17 = self.geom.add_line_loop([fsp_l139, -fsp_l118, fsp_l117, fsp_l120, -fsp_l137])
        fsp_ll18 = self.geom.add_line_loop([fsp_l119, fsp_l123, -fsp_l142, -fsp_l140, -fsp_l120])
        fsp_ll19 = self.geom.add_line_loop([fsp_l121, fsp_l123, -fsp_l189, -fsp_l161, fsp_l160, fsp_l165, fsp_l193, -fsp_l122])
        fsp_ll20 = self.geom.add_line_loop([fsp_l189, -fsp_l142, -fsp_l134, fsp_l135, fsp_l182])
        fsp_ll21 = self.geom.add_line_loop([fsp_l161, -fsp_l182, -fsp_l181, -fsp_l157, fsp_l156])
        fsp_ll22 = self.geom.add_line_loop([fsp_l99, fsp_l102, -fsp_l133, fsp_l135, -fsp_l181, -fsp_l100])
        fsp_ll23 = self.geom.add_line_loop([fsp_l101, -fsp_l103, -fsp_l85, -fsp_l81, fsp_l83, fsp_l97, fsp_l99])
        fsp_ll24 = self.geom.add_line_loop([fsp_l97, fsp_l100, -fsp_l157, -fsp_l154, -fsp_l98])
        fsp_ll25 = self.geom.add_line_loop([fsp_l82, fsp_l87, fsp_l125, -fsp_l153, -fsp_l98, -fsp_l83])
        fsp_ll26 = self.geom.add_line_loop([fsp_l154, fsp_l156, fsp_l160, fsp_l164, -fsp_l162, -fsp_l58, -fsp_l155, -fsp_l153])
        fsp_ll27 = self.geom.add_line_loop([fsp_l124, fsp_l127, fsp_l186, -fsp_l159, -fsp_l155, -fsp_l125])
        fsp_ll28 = self.geom.add_line_loop([fsp_l81, fsp_l84, -fsp_l86, -fsp_l82])
        fsp_ll29 = self.geom.add_line_loop([fsp_l85, fsp_l104, -fsp_l198, -fsp_l131, -fsp_l88, -fsp_l84])
        fsp_ll30 = self.geom.add_line_loop([fsp_l105, -fsp_l177, fsp_l176, fsp_l180, fsp_l200, -fsp_l104, fsp_l103])
        fsp_ll31 = self.geom.add_line_loop([fsp_l158, fsp_l163, -fsp_l191, -fsp_l187, -fsp_l159])
        fsp_ll32 = self.geom.add_line_loop([fsp_l151, -fsp_l187, -fsp_l186, -fsp_l149, fsp_l148])
        fsp_ll33 = self.geom.add_line_loop([fsp_l149, -fsp_l127, fsp_l126, fsp_l129, -fsp_l147])
        fsp_ll34 = self.geom.add_line_loop([fsp_l144, -fsp_l150, -fsp_l147, fsp_l148, -fsp_l146, -fsp_l143])
        fsp_ll35 = self.geom.add_line_loop([fsp_l124, fsp_l126, fsp_l128, -fsp_l130, -fsp_l88, -fsp_l86, fsp_l87])
        fsp_ll36 = self.geom.add_line_loop([fsp_l130, fsp_l132, -fsp_l190, -fsp_l175, fsp_l174, fsp_l179, fsp_l196, -fsp_l131])
        fsp_ll37 = self.geom.add_line_loop([fsp_l128, fsp_l132, -fsp_l152, -fsp_l150, -fsp_l129])
        fsp_ll38 = self.geom.add_line_loop([fsp_l198, -fsp_l200, -fsp_l197, fsp_l196])
        fsp_ll39 = self.geom.add_line_loop([fsp_l162, fsp_l166, fsp_l199, -fsp_l113, fsp_l112, fsp_l114, -fsp_l163])
        fsp_ll40 = self.geom.add_line_loop([fsp_l178, fsp_l180, -fsp_l197, -fsp_l179])
        fsp_ll41 = self.geom.add_line_loop([fsp_l164, fsp_l166, -fsp_l194, -fsp_l165])
        fsp_ll42 = self.geom.add_line_loop([fsp_l195, -fsp_l199, -fsp_l194, fsp_l193])
        
        fsp_s1 = self.geom.add_plane_surface(fsp_ll1, holes=None)
        fsp_s2 = self.geom.add_plane_surface(fsp_ll2, holes=None)
        fsp_s3 = self.geom.add_plane_surface(fsp_ll3, holes=None)
        fsp_s4 = self.geom.add_plane_surface(fsp_ll4, holes=None)
        fsp_s5 = self.geom.add_plane_surface(fsp_ll5, holes=None)
        fsp_s6 = self.geom.add_plane_surface(fsp_ll6, holes=None)
        fsp_s7 = self.geom.add_plane_surface(fsp_ll7, holes=None)
        fsp_s8 = self.geom.add_plane_surface(fsp_ll8, holes=None)
        fsp_s9 = self.geom.add_plane_surface(fsp_ll9, holes=None)
        fsp_s10 = self.geom.add_plane_surface(fsp_ll10, holes=None)
        fsp_s11 = self.geom.add_plane_surface(fsp_ll11, holes=None)
        fsp_s12 = self.geom.add_plane_surface(fsp_ll12, holes=None)
        fsp_s13 = self.geom.add_plane_surface(fsp_ll13, holes=None)
        fsp_s14 = self.geom.add_plane_surface(fsp_ll14, holes=None)
        fsp_s15 = self.geom.add_plane_surface(fsp_ll15, holes=None)
        fsp_s16 = self.geom.add_plane_surface(fsp_ll16, holes=None)
        fsp_s17 = self.geom.add_plane_surface(fsp_ll17, holes=None)
        fsp_s18 = self.geom.add_plane_surface(fsp_ll18, holes=None)
        fsp_s19 = self.geom.add_plane_surface(fsp_ll19, holes=None)
        fsp_s20 = self.geom.add_plane_surface(fsp_ll20, holes=None)
        fsp_s21 = self.geom.add_plane_surface(fsp_ll21, holes=None)
        fsp_s22 = self.geom.add_plane_surface(fsp_ll22, holes=None)
        fsp_s23 = self.geom.add_plane_surface(fsp_ll23, holes=None)
        fsp_s24 = self.geom.add_plane_surface(fsp_ll24, holes=None)
        fsp_s25 = self.geom.add_plane_surface(fsp_ll25, holes=None)
        fsp_s26 = self.geom.add_plane_surface(fsp_ll26, holes=None)
        fsp_s27 = self.geom.add_plane_surface(fsp_ll27, holes=None)
        fsp_s28 = self.geom.add_plane_surface(fsp_ll28, holes=None)
        fsp_s29 = self.geom.add_plane_surface(fsp_ll29, holes=None)
        fsp_s30 = self.geom.add_plane_surface(fsp_ll30, holes=None)
        fsp_s31 = self.geom.add_plane_surface(fsp_ll31, holes=None)
        fsp_s32 = self.geom.add_plane_surface(fsp_ll32, holes=None)
        fsp_s33 = self.geom.add_plane_surface(fsp_ll33, holes=None)
        fsp_s34 = self.geom.add_plane_surface(fsp_ll34, holes=None)
        fsp_s35 = self.geom.add_plane_surface(fsp_ll35, holes=None)
        fsp_s36 = self.geom.add_plane_surface(fsp_ll36, holes=None)
        fsp_s37 = self.geom.add_plane_surface(fsp_ll37, holes=None)
        fsp_s38 = self.geom.add_plane_surface(fsp_ll38, holes=None)
        fsp_s39 = self.geom.add_plane_surface(fsp_ll39, holes=None)
        fsp_s40 = self.geom.add_plane_surface(fsp_ll40, holes=None)
        fsp_s41 = self.geom.add_plane_surface(fsp_ll41, holes=None)
        fsp_s42 = self.geom.add_plane_surface(fsp_ll42, holes=None)

        
        surfaces = [fsp_s1, fsp_s2, fsp_s3, fsp_s4, fsp_s5, fsp_s6, fsp_s7, fsp_s8, fsp_s9, fsp_s10, fsp_s11, fsp_s12, fsp_s13, fsp_s14, fsp_s15, fsp_s16, fsp_s17, fsp_s18, fsp_s19, fsp_s20, fsp_s21, fsp_s22, fsp_s23, fsp_s24, fsp_s25, fsp_s26, fsp_s27, fsp_s28, fsp_s29, fsp_s30, fsp_s31, fsp_s32, fsp_s33, fsp_s34, fsp_s35, fsp_s36, fsp_s37, fsp_s38, fsp_s39, fsp_s40, fsp_s41, fsp_s42]
        
        return surfaces


    def surfaceloop(self):
        return self.geom.add_surface_loop(self.surfaces())
        
    def phys_surface(self):        
        return self.geom.add_physical(self.surfaces())
        
    def volume(self):
        return self.geom.add_volume(self.surfaceloop())
        
    def phys_volume(self):
        return self.geom.add_physical(self.volume())
  
#mesh = pygmsh.generate_mesh(geom)

#meshio.write("spherical.vtk", mesh)


  
#mesh = pygmsh.generate_mesh(geom)

#meshio.write("spherical.vtk", mesh)





