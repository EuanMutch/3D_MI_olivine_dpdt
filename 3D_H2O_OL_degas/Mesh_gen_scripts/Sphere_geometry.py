
import math as m
import numpy as np


# functions for calculating plane-sphere intersection


def pl_eq(p1, p2, p3):
    v1 = p2-p1
    v2 = p3-p1

    v3 = np.cross(v1, v2)

    A, B, C = v3[0], v3[1], v3[2]

    D = -A*p1[0] - B*p1[1] - C*p1[2]

    return A, B, C, D


def pl_sp_int(A, B, C, D, xs, ys, zs, r):

    xc = xs - ((A*(A*xs + B*ys + C*zs + D))/(A**2 + B**2 + C**2))
    yc = ys - ((B*(A*xs + B*ys + C*zs + D))/(A**2 + B**2 + C**2))
    zc = zs - ((C*(A*xs + B*ys + C*zs + D))/(A**2 + B**2 + C**2))

    d = m.sqrt((xs - xc)**2 + (ys - yc)**2 + (zs - zc)**2)

    # Conditions of intersection

    if r > d:
        return 2 #Plane intersection

    elif r == d:
        return 1 #Plane tangent

    else:
        return 0 # No intersection

    
def pl_pos(A, B, C, D, xs, ys, zs): 
    dd = A*xs + B*ys + C*zs + D
    return dd    
