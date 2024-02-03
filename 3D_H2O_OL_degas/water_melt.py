



# Python script to calculate H2O diffusivities in melt for different compositions. 

import pandas as pd
import numpy as np
import math as m
from dolfin import *


# Import glass data

#df1 = pd.read_csv('Seguam_melt.csv')

# Calculate mole fractions

def Dwtrdf(df, T, P, spec):

    # Input parameters
    # df - pandas dataframe of oxides corresponding to melt composition
    # T - temperature (C) - float
    # P - pressure (GPa) - float
    # spec - diffusion coefficient of required species - str - H2Om for molecular water, OH for hydroxyl and H2Ot for total water.

    oxides = np.array(['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'MgO', 'MnO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O'])

    oxide_masses = {'SiO2': 28.085 + (2*15.999) , 'TiO2': 47.867 + (2*15.999), 'Al2O3': (2*26.982) + (3*15.999), 'Fe2O3': (2*55.845) + (3*15.999), 'FeO': 55.845 + 15.999 , 'MgO': 24.305 + 15.999, 'MnO': 54.938 + 15.999, 'CaO': 40.078 + 15.999, 'Na2O': (2*22.99) + 15.999, 'K2O': (2*39.098) + 15.999, 'P2O5': (5*30.974) + (5*15.999), 'Cr2O3': (2*51.996) + (3*15.999), 'NiO': 58.693 + 15.999, 'CoO': 59.933 + 15.999, 'H2O': (2*1.008) + 15.999} 

    mf = {} 
    
    # Calculate moles of each element - easier to do by columns

    for zz in oxides:
        mf['M_{0}'.format(zz)] = df['{0}'.format(zz)].values[0]/oxide_masses['{0}'.format(zz)]

    mf['M_AlO1.5'] = mf['M_Al2O3']*2
    mf['M_FeO1.5'] = mf['M_Fe2O3']*2
    mf['M_CrO1.5'] = mf['M_Cr2O3']*2
    mf['M_NaO0.5'] = mf['M_Na2O']*2
    mf['M_KO0.5'] = mf['M_K2O']*2
    mf['M_PO2.5'] = mf['M_P2O5']*2
    mf['M_HO0.5'] = mf['M_H2O']*2

    # Calculate the sum of moles

    m_tot = mf['M_SiO2'] + mf['M_TiO2'] + mf['M_AlO1.5'] + mf['M_CrO1.5'] + mf['M_FeO'] + mf['M_FeO1.5'] + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_NaO0.5'] + mf['M_KO0.5'] + mf['M_PO2.5']

    m_totO = (2*mf['M_SiO2']) + (2*mf['M_TiO2']) + (3*mf['M_Al2O3']) + (3*mf['M_Cr2O3']) + mf['M_FeO'] + (3*mf['M_Fe2O3']) + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_Na2O'] + mf['M_K2O'] + (5*mf['M_P2O5'])

    T_K = T + 273

    # Molar mass of a single oxygen basis g/mol
    W = df['Total'].values[0]/m_totO 

    # Mole fraction of Si among all cations
    Y_Si = mf['M_SiO2']/m_tot 

    # Mole fraction of H2Ot on a single oxygen basis
    X = (mf['M_H2O'])/((mf['M_H2O']) + ((100.0-df['H2O'].values[0])/W)) 

    #Equilibrium constant of speciation reaction    
    K = m.exp(2.6*Y_Si - 4339.0*Y_Si/T_K) 

    dXOH_2dX = (0.5-X)/m.sqrt(X*(1-X)*(4.0/K-1.0)+0.25)

    dXm_dX = 1.0 - dXOH_2dX

    alpha = -94.07 + 74.112*Y_Si + (198508.0-166674.0*Y_Si)/T_K

    DOH_D0 = m.exp(-56.09-115.93*Y_Si + 160.54*Y_Si**0.5 - 3970.0*Y_Si**0.5/T_K)

    #pre-exponential factor
    D0 = m.exp(8.02-31*Y_Si+2.348*P*Y_Si+(121824*Y_Si-118323*Y_Si**0.5-(10016*Y_Si-3648)*P)/T_K) 

    #diffusivity of molecular water
    D_H2O = D0*m.exp(alpha*X) 

    #diffusivity of hydroxyl group
    D_OH = DOH_D0*D0 

    #diffusivity of total water
    D_H2Ot = D_H2O*dXm_dX + D_OH*dXOH_2dX 

    if spec == "H2Om":
        return D_H2O

    elif spec == "OH":
        return D_OH

    elif spec == "H2Ot":
        return D_H2Ot

#################################################################

# Version of function with input for H2O (not from imported dataframe)

def Dwtri(C, df, T, P, spec):

    # Input parameters
    # C - water concentration (in wt%)
    # df - pandas dataframe of oxides corresponding to melt composition
    # T - temperature (C) - float
    # P - pressure (GPa) - float
    # spec - diffusion coefficient of required species - str - H2Om for molecular water, OH for hydroxyl and H2Ot for total water.

    oxides = np.array(['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'MgO', 'MnO', 'CaO', 'Na2O', 'K2O', 'P2O5'])

    oxide_masses = {'SiO2': 28.085 + (2*15.999) , 'TiO2': 47.867 + (2*15.999), 'Al2O3': (2*26.982) + (3*15.999), 'Fe2O3': (2*55.845) + (3*15.999), 'FeO': 55.845 + 15.999 , 'MgO': 24.305 + 15.999, 'MnO': 54.938 + 15.999, 'CaO': 40.078 + 15.999, 'Na2O': (2*22.99) + 15.999, 'K2O': (2*39.098) + 15.999, 'P2O5': (5*30.974) + (5*15.999), 'Cr2O3': (2*51.996) + (3*15.999), 'NiO': 58.693 + 15.999, 'CoO': 59.933 + 15.999, 'H2O': (2*1.008) + 15.999} 

    mf = {} 
    
    # Calculate moles of each element - easier to do by columns

    for zz in oxides:
        mf['M_{0}'.format(zz)] = df['{0}'.format(zz)].values[0]/oxide_masses['{0}'.format(zz)]

    mf['M_H2O'] = C/oxide_masses['H2O']
    mf['M_AlO1.5'] = mf['M_Al2O3']*2
    mf['M_FeO1.5'] = mf['M_Fe2O3']*2
    mf['M_CrO1.5'] = mf['M_Cr2O3']*2
    mf['M_NaO0.5'] = mf['M_Na2O']*2
    mf['M_KO0.5'] = mf['M_K2O']*2
    mf['M_PO2.5'] = mf['M_P2O5']*2
    mf['M_HO0.5'] = mf['M_H2O']*2

    # Calculate the sum of moles

    m_tot = mf['M_SiO2'] + mf['M_TiO2'] + mf['M_AlO1.5'] + mf['M_CrO1.5'] + mf['M_FeO'] + mf['M_FeO1.5'] + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_NaO0.5'] + mf['M_KO0.5'] + mf['M_PO2.5']

    m_totO = (2*mf['M_SiO2']) + (2*mf['M_TiO2']) + (3*mf['M_Al2O3']) + (3*mf['M_Cr2O3']) + mf['M_FeO'] + (3*mf['M_Fe2O3']) + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_Na2O'] + mf['M_K2O'] + (5*mf['M_P2O5'])

    T_K = T + 273

    # Molar mass of a single oxygen basis g/mol
    W = df['Total'].values[0]/m_totO 

    # Mole fraction of Si among all cations
    Y_Si = mf['M_SiO2']/m_tot 

    # Mole fraction of H2Ot on a single oxygen basis
    X = (mf['M_H2O'])/((mf['M_H2O']) + ((100.0-C)/W)) 

    #Equilibrium constant of speciation reaction    
    K = m.exp(2.6*Y_Si - 4339.0*Y_Si/T_K) 

    dXOH_2dX = (0.5-X)/m.sqrt(X*(1-X)*(4.0/K-1.0)+0.25)

    dXm_dX = 1.0 - dXOH_2dX

    alpha = -94.07 + 74.112*Y_Si + (198508.0-166674.0*Y_Si)/T_K

    DOH_D0 = m.exp(-56.09-115.93*Y_Si + 160.54*Y_Si**0.5 - 3970.0*Y_Si**0.5/T_K)

    #pre-exponential factor
    D0 = m.exp(8.02-31*Y_Si+2.348*P*Y_Si+(121824*Y_Si-118323*Y_Si**0.5-(10016*Y_Si-3648)*P)/T_K) 

    #diffusivity of molecular water
    D_H2O = D0*m.exp(alpha*X) 

    #diffusivity of hydroxyl group
    D_OH = DOH_D0*D0 

    #diffusivity of total water
    D_H2Ot = D_H2O*dXm_dX + D_OH*dXOH_2dX 

    if spec == "H2Om":
        return D_H2O

    elif spec == "OH":
        return D_OH

    elif spec == "H2Ot":
        return D_H2Ot

#################################################################

# Version of function with input for H2O of equilibrium olivine composition

def Dwtrol(Col, Kd, df, T, P, spec):

    # Input parameters
    # Col - water concentration of olivine in equilibrium with the melt inclusion (in ppm)
    # Kd - olivine/melt partition coefficient
    # df - pandas dataframe of oxides corresponding to melt composition of inclusion
    # T - temperature (C) - float
    # P - pressure (GPa) - float
    # spec - diffusion coefficient of required species - str - H2Om for molecular water, OH for hydroxyl and H2Ot for total water.

    # Calculate composition of melt in equilibrium with olivine and convert to wt%
    Cmelt = (Col/Kd)/10000.0 

    oxides = np.array(['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'MgO', 'MnO', 'CaO', 'Na2O', 'K2O', 'P2O5'])

    oxide_masses = {'SiO2': 28.085 + (2*15.999) , 'TiO2': 47.867 + (2*15.999), 'Al2O3': (2*26.982) + (3*15.999), 'Fe2O3': (2*55.845) + (3*15.999), 'FeO': 55.845 + 15.999 , 'MgO': 24.305 + 15.999, 'MnO': 54.938 + 15.999, 'CaO': 40.078 + 15.999, 'Na2O': (2*22.99) + 15.999, 'K2O': (2*39.098) + 15.999, 'P2O5': (5*30.974) + (5*15.999), 'Cr2O3': (2*51.996) + (3*15.999), 'NiO': 58.693 + 15.999, 'CoO': 59.933 + 15.999, 'H2O': (2*1.008) + 15.999} 

    mf = {} 
    
    # Calculate moles of each element - easier to do by columns

    for zz in oxides:
        mf['M_{0}'.format(zz)] = df['{0}'.format(zz)].values[0]/oxide_masses['{0}'.format(zz)]

    mf['M_H2O'] = Cmelt/oxide_masses['H2O']
    mf['M_AlO1.5'] = mf['M_Al2O3']*2
    mf['M_FeO1.5'] = mf['M_Fe2O3']*2
    mf['M_CrO1.5'] = mf['M_Cr2O3']*2
    mf['M_NaO0.5'] = mf['M_Na2O']*2
    mf['M_KO0.5'] = mf['M_K2O']*2
    mf['M_PO2.5'] = mf['M_P2O5']*2
    mf['M_HO0.5'] = mf['M_H2O']*2

    # Calculate the sum of moles

    m_tot = mf['M_SiO2'] + mf['M_TiO2'] + mf['M_AlO1.5'] + mf['M_CrO1.5'] + mf['M_FeO'] + mf['M_FeO1.5'] + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_NaO0.5'] + mf['M_KO0.5'] + mf['M_PO2.5']

    m_totO = (2*mf['M_SiO2']) + (2*mf['M_TiO2']) + (3*mf['M_Al2O3']) + (3*mf['M_Cr2O3']) + mf['M_FeO'] + (3*mf['M_Fe2O3']) + mf['M_MgO'] + mf['M_MnO'] + mf['M_CaO'] + mf['M_Na2O'] + mf['M_K2O'] + (5*mf['M_P2O5'])

    T_K = T + 273

    # Molar mass of a single oxygen basis g/mol
    W = df['Total'].values[0]/m_totO 

    # Mole fraction of Si among all cations
    Y_Si = mf['M_SiO2']/m_tot 

    # Mole fraction of H2Ot on a single oxygen basis
    X = (mf['M_H2O'])/((mf['M_H2O']) + ((100.0-Cmelt)/W)) 

    #Equilibrium constant of speciation reaction    
    K = m.exp(2.6*Y_Si - 4339.0*Y_Si/T_K) 

    dXOH_2dX = (0.5-X)/sqrt(X*(1-X)*(4.0/K-1.0)+0.25)

    dXm_dX = 1.0 - dXOH_2dX

    alpha = -94.07 + 74.112*Y_Si + (198508.0-166674.0*Y_Si)/T_K

    DOH_D0 = m.exp(-56.09-115.93*Y_Si + 160.54*Y_Si**0.5 - 3970.0*Y_Si**0.5/T_K)

    #pre-exponential factor
    D0 = m.exp(8.02-31*Y_Si+2.348*P*Y_Si+(121824*Y_Si-118323*Y_Si**0.5-(10016*Y_Si-3648)*P)/T_K) 

    #diffusivity of molecular water
    D_H2O = D0*exp(alpha*X) 

    #diffusivity of hydroxyl group
    D_OH = DOH_D0*D0 

    #diffusivity of total water
    D_H2Ot = D_H2O*dXm_dX + D_OH*dXOH_2dX 

    if spec == "H2Om":
        return D_H2O

    elif spec == "OH":
        return D_OH

    elif spec == "H2Ot":
        return D_H2Ot


