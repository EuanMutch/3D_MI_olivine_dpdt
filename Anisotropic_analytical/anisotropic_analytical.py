
# Script to calculate analytical solution for the anisotropic diffusive equilibration of melt inclusion along a degassing pathway. Based on equations 14-16 in the main manuscript. 


# Import python modules
import matplotlib.pyplot as plt
import math as m
import pandas as pd
import numpy as np
from scipy.special import elliprf
from scipy.integrate import cumtrapz

import os
import sys
import time
import datetime


def Aniso_analytical(mi_r, f_Kd, f_aniso, f_dpdt, f_degas, p_degas):     
    """
    Function for calculating analytical solution for anisotropic equilibration of a melt inclusion. Returns the composition of the melt inclusion following ascent upon degassing.

    ...

    Parameters
    ----------
    mi_r : float
        Melt inclusion radius (microns)
    f_Kd : float
        olivine-melt partition coefficient
    f_aniso : float
        Diffusive anisotropy in the olivine
    f_dpdt : float
        Decompression rate (MPa s-1)
    f_degas : array
        Water content along the degassing pathway (wt%)
    p_degas : array
        Pressure along the degassing pathway (MPa)

    Returns
    -------
    cmi : array
        Water contents in the centre of the melt inclusion (wt%)
    """
    
    c = f_degas
    p = p_degas

    t = (p[0] - p) / f_dpdt
    
    # DIFFUSION COEFFICIENTS

    # Olivine diffusion coefficient - H+ with anisotropy - values from Barth et al. (2019)
    D_100 = ((9.6e-6)*m.exp(-125000.0/(8.314*T_K)))*(1e12)
    D_010 = (1.0/f_aniso)*D_100 # Assume only 10x anisotropy for now. 40 
    D_001 = (1.0/f_aniso)*D_100

    dr_olm = 1.2 # density ratio of olivine/melt

    # Set the partition coefficient
    Kd = f_Kd

    D_eff = np.sqrt(D_100 * D_010 * D_001) / elliprf(D_100**-1.0, D_010**-1.0, D_001**-1.0)

    tau = (1.0 / (dr_olm * f_Kd)) * mi_r**2 / (3.0 * D_eff)
      
    cmi = np.exp(-t / tau) * (cumtrapz(c * np.exp(t / tau) / tau, x=t, initial=0.0) + c[0]) 
    
    return cmi
    
       
# Import files 


dfm = pd.read_csv('Seguam_degas.csv')

h2o_degas = dfm['H2O'].values
P_degas = dfm['P_MPa'].values

# define parameters

kd = 0.000459 # Partition coefficient between olivine and melt
aniso = 15.6 # Anisotropy in olivine
T_K = 1070.0 + 273.0 # Temperature in K
mi_r = 30.0 # Melt inclusion radius
dpdt = 0.02 # decompression rate MPa/s

C_end = Aniso_analytical(mi_r, kd, aniso, dpdt, h2o_degas, P_degas)

plt.plot(P_degas, h2o_degas)
plt.plot(P_degas, C_end)

plt.show()



