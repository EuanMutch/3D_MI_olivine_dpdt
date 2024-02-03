#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 16:01:46 2020

@author: ejfm2
"""

#import numpy as np
import cython

from libc.math cimport exp
from libc.math cimport log
from libc.math cimport pow
from libc.math cimport sqrt
from libc.math cimport copysignf



cdef double FNA(double TK):
    cdef double x = (166800000.0 - 193080.0 * (TK - 273.15) + 186.4 * pow(TK - 273.15, 2.0) - 0.071288 * (pow((TK - 273.15), 3.0))) * 1.01325
    return x

cdef double FNB(double TK):
    cdef double x = 1.01325 * (73030000.0 - 71400.0 * (TK - 273.15) + 21.57 * pow((TK - 273.15), 2.0))
    return x

@cython.cdivision(True)    
cdef double FNC(double TK):
    cdef double R = 83.14321
    cdef double x = 1.01325 * (exp(-11.071 + 5953.0 / TK - 2746000.0 / pow(TK, 2.0) + 464600000.0 / pow(TK, 3.0)) * 0.5 * R * R * pow(TK, 2.5) / 1.02668 + 40123800.0)
    return x

@cython.cdivision(True)     
cdef double FNF(double V, double TK, double A, double B, double P):
    cdef double R = 83.14321
    cdef double x = R * TK / (V - B) - A / ((V * V + B * V) * sqrt(TK)) - P
    return x

@cython.cdivision(True) 
cdef (double, double) MRK(double P, double TK): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
    cdef double R = 83.14321
    cdef double B_1 = 14.6
    cdef double B_2 = 29.7
    cdef int X_1
    cdef double B
    cdef double A
    cdef double Temp2
    cdef double Q
    cdef double Temp1
    cdef double F_1
    cdef double F_2
    cdef double V
    cdef double G_1
    cdef double G_2
    cdef double fCO2o
    cdef double fH2Oo
    cdef (double, double) output 
    for X_1 in range(2): #loops twice, once for each CO2 and H2O
        B = X_1 * B_1 + (1 - X_1) * B_2
        A = X_1**2.0 * FNA(TK) + 2.0 * X_1 * (1.0 - X_1) * FNC(TK) + (1.0 - X_1)**2.0 * FNB(TK)
        Temp2 = B + 5.0
        Q = 1.0
        Temp1 = 0.0
        while abs(Temp2 - Temp1) >= 0.00001:
            Temp1 = Temp2
            F_1 = (FNF(Temp1 + 0.01, TK, A, B, P) - FNF(Temp1, TK, A, B, P)) / 0.01
            Temp2 = Temp1 - Q * FNF(Temp1, TK, A, B, P) / F_1
            F_2 = (FNF(Temp2 + 0.01, TK, A, B, P) - FNF(Temp2, TK, A, B, P)) / 0.01
            if F_2 * F_1 <= 0:
                Q = Q / 2.
            if abs(Temp2 - Temp1) > 0.00001:
                F_1 = F_2
        V = Temp2
        G_1 = log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * FNA(TK) + (1 - X_1) * FNC(TK)) * log((V + B) / V) / (R * TK**1.5 * B)
        G_1 = G_1 + (log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - log(P * V / (R * TK))
        G_1 = exp(G_1)
        G_2 = log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * FNC(TK) + (1 - X_1) * FNB(TK)) * log((V + B) / V) / (R * TK**1.5 * B)
        G_2 = G_2 + (log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - log(P * V / (R * TK))
        G_2 = exp(G_2)
        if X_1 == 0:
            fCO2o = G_2 * P #The fugacity of CO2
        if X_1 == 1:
            fH2Oo = G_1 * P #The fugacity of H2O
    output = fCO2o, fH2Oo
    return output

@cython.cdivision(True) 
cdef list SatPress(int M, int comp, double SiO2, double WtH2O, double PPMCO2, double TK, double P):
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Input variables
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #routine is a list, the first element is M
    #M = 0 (calc saturation pressure), 1 (equilibrium speciation), 2 (isobar calculation)
    #if M = 1 or 2, the second element in M is a starting pressure
    #comp = 0 (basalt), 1 (rhyolite)
    #SiO2 = SiO2 content in weight percent
    #WtH2O = H2O content in weight percent
    #PPMCO2 = CO2 content in ppm
    #TK = temperature in kelvin
    cdef double press
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #local variables
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #If vapor saturation pressure calculation, routine[0] = 1
    if M != 0: #set pressure cases other than vapor sat pressure calculation
        press = P
    
    cdef double Z = 10. #This is the increment of pressure, default is 10 bar increment
    
    cdef double H2Ocoeff = (-0.0000304 + 0.00000129 * SiO2) / 0.0000328
    cdef double CO2coeff = (0.0000087 - 0.0000001698 * SiO2) / 0.00000038
    cdef double temp
    cdef double Y
    cdef double GH2O
    cdef double GCO2
    cdef double xb
    cdef double xH2Om
    cdef double SNCO2
    cdef double SNH2O
    cdef double SNO
    cdef double xCO2m
    cdef double XOHm
    cdef double XOHmo
    cdef double XH2O
    cdef double XO
    cdef double derath
    cdef double fx
    cdef double fxp
    cdef double fCO2o
    cdef double fH2Oo
    cdef double Yo
    cdef double changer
    
    
    if M == 0: #if M == 1, pressure is approximately known
        temp = WtH2O + PPMCO2 / 250.
        if temp < 0.5:
            press = 20.
            Z = 1.
        elif temp < 1:
            press = 40.
            Z = 20.
        elif temp < 2:
            press = 100.
            Z = 100.
        elif temp < 3:
            press = 500.
            Z = 100.
        elif temp < 4:
            press = 1000.
            Z = 100.
        elif temp < 5:
            press = 1500.
            Z = 100.
        elif temp < 6:
            press = 2000.
            Z = 100.
        else:
            press = 3000.
            Z = 100.
    if M == 2: #M == 2 is for isobars, pressure is not variable
        changer = WtH2O #H2O is the variable
        Z = 0.1
    else:
        changer = press #pressure is the variable
    
    #initialize Y for Newton's method
    Y = 1.0
    GH2O = 2.0
    GCO2 = 0.0
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #isobar loop
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    '''
    Uses Newton's method to move changer (WtH2O or P, depending on calling subroutine, and thus M).  Looping stops when
    the two mole fractions add up to 1.
    '''
    while abs(GH2O + GCO2 - 1) >= 0.0001:
        #For low-H2O samples, XH2Om is a simple function.  Avoids problems in code resulting in division by zero.
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Calculate xH2Om,xCO2m
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if WtH2O < 0.6:
            xb = 0.0
            XOHm = 0.0
            if comp == 0: #if the composition is basaltic
                xH2Om = exp(-5.827) * WtH2O**1.855
                SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                SNH2O = WtH2O / 18.015 #relative moles of H2O
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O) #mol fraction of CO2 relative to H2O+CO2+O
                
            else: #if the composition is rhyolitic
                xH2Om = exp(-5.526) * WtH2O**1.977
                SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                SNH2O = WtH2O / 18.015 #relative moles of H2O
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O) #mol fraction of CO2 relative to H2O+CO2+O
        
        #For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
        else: #For WtH2O > 0.6, use method of Silver, Newman, Stolper et al.
            XOHm = 0.01
            XOHmo = XOHm
            SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
            SNH2O = WtH2O / 18.015 #relative moles of H2O
            if comp == 0: #if composition is basaltic
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
            else: #if composition is rhyolitic
                SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
            XH2O = SNH2O / (SNH2O + SNO) #mol fraction of water, relative to H2O+O
            XO = 1 - XH2O #mol fraction of O, relative to H2O+O
            xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)
            if comp == 0: #if composition is basaltic
                xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
            else: #if composition is rhyolitic
                xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
            #The following three if/then statements prevent infinite loops and division by zero
            if XH2O - 0.6 * XOHm <= 0:
                XOHm = (0.999 * 2 * XH2O + XOHmo) / 2
            if XO - 0.5 * XOHm <= 0:
                XOHm = (0.999 * 2 * XO + XOHmo) / 2
            if XOHm**2 == 0:
                XOHm = XOHmo / 2
            derath = 1.0
            if comp == 0: #if composition is basaltic
                while abs(derath) >= 0.00001: #loop to define XOHm
                    fx = 9.143 - 3.295 * (XOHm - 1) - 2 * 6.019 * (XO - XOHm) - 2 * 0.572 * (XH2O - XOHm) + log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                    fxp = -3.295 + 2 * (6.019 + 0.572) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * ((XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                    derath = fx / fxp
                    XOHm = XOHm - derath
            else: #if composition is rhyolitic
                while abs(derath)>= 0.00001: #loop to define XOHm
                    fx = 9.345 - 4.304 * (XOHm - 1) - 2 * 6.277 * (XO - XOHm) - 2 * 2.328 * (XH2O - XOHm) + log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                    fxp = -4.304 + 2 * (6.277 + 2.328) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * ((XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                    derath = fx / fxp
                    XOHm = XOHm - derath
            xH2Om = xb - 0.5 * XOHm
            
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Calculate gas
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        fCO2o,fH2Oo = MRK(press,TK) #Calls MRK procedure to define the endmember fugacities of H2O and CO2
        '''
        GH2O and GCO2 are the mol fractions of H2O and CO2 in the gas phase.  They are dependent on XCO2m and XH2Om
        as well as T and P and endmember H2O and CO2 fugacities.
        Standard state for basalt is 1 bar and 1200C.  All values, plus Delta V and Delta H from Dixon et al. (1995).
        Solubilities at standard state are 0.5 ppm CO2 and 0.11 wt.% H2O for a basalt with 49% SiO2. These values
        are equivalent to mol fractions of .00000038 and .0000328, respectively.
        Solubilities vary accoding to wt.% SiO2 in accord with H2Ocoeff and CO2coeff derived by Dixon (1997) and
        slightly modified so that a basalt with 49% SiO2 is equivalent to the tholeiite in the original (1995) model.
        
        If GH2O+GCO2 sum to > 1, then pressure must be decreased or either Wt.H2O or PPMCO2 need to be decreased.
        '''
        
        if comp == 0: #if composition is basaltic
            GH2O = (xH2Om / 0.0000328 / H2Ocoeff) * (1 / fH2Oo) / exp(-12 * (press - 1) / (41.84 * 1.9872 * TK))
            GCO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / exp(-23 * (press - 1) / (41.84 * 1.9872 * TK))
        #
        #Standard state for water in rhyolite is 799 bars and 850C.  Solubility under those conditions is 3.45 wt.% _
        #which corresponds to xH2Om of 0.0323.  These values and partial molar volume of 10 cc/mol from Silver (1977).
        #Delta H of -4420 kcal/mol is for H2O dissolution in albite from Silver et al. (1990).
        #
        #Standard state for CO2 in rhyolite is 750 bars and 850C.  Solubility under those conditions is 340 ppm _
        #which corresponds to xCO2m of 0.000397.  Partial molar volume is 28 cc/mol. Values from Blank et al. (1993).
        #DeltaH for CO2 dissolution in melt (4861 kcal/mol) from Fogel and Rutherford (1990).
        #
        else: #if composition is rhyolitic
            GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / exp(-5 * (press - 799) / (41.84 * 1.9872 * TK) + 4420 / 1.9872 * (1 / TK - 1 / 1123.16))
            GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / exp(-28 * (press - 750) / (41.84 * 1.9872 * TK) + 4861 * (1 / TK - 1 / 1123.16) / 1.9872)
        
        #Redefine variables
        if abs(GH2O + GCO2 - 1) >= 0.0001:
            Yo = Y
            Y = copysignf(1.0, (-1. * (GH2O + GCO2 - 1)))
            if Y + Yo == 0:
                Z = Z / 2.
            if M == 2: #Isobaric: variable is H2O
                changer = changer + Z * Y / 2.
                WtH2O = changer
            if M != 2: #Isocompositional: variable is pressure
                changer = changer - Z * Y / 2.
                press = changer
    
    if comp == 0: #if composition is basaltic
        WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% molecular H2O in melt
        WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 36.594) #Calculate Wt% hydroxyl in melt
    else: #if composition is rhyolitic
        WtH2Om = (xH2Om * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% molecular H2O in melt
        WtOHm = (0.5 * XOHm * 1801.5) / ((xb * 18.015) + (1 - xb) * 32.5) #Calculate Wt% hydroxyl in melt
        
    if WtH2O < 0.6:
        WtOHm = WtH2O - WtH2Om
    
    cdef list output = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
    return output




###########################################################################
# Degassing path
###########################################################################



@cython.cdivision(True) 
cpdef degas(int comp, double WtH2O, double PPMCO2, double SiO2, double T, int style, double excess, int steps):

    cdef double TK = T + 273.0
    cdef int routine = 0 #routine = [0] means that SatPress function varies pressure
    cdef list out
    cdef list output
    cdef double press
    cdef double GH2O
    cdef double GCO2
    cdef double WtH2Olast
    cdef double inc
    cdef double Drop
    cdef int seq
    cdef int loopcount
    cdef double CarbonVFrac
    cdef double XX
    cdef double PPMCO2o
    cdef double WTH2Oo
    cdef double CarbonVapor
    cdef double WaterVapor
    cdef double Xz
    cdef double WFCO2v
    cdef double WtVap
    cdef double MassVCO2
    cdef double MassVH2O
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Test input variables for errors
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#==============================================================================
#     if comp == 0:
#         if SiO2 < 40 or SiO2 > 49:
#             error[0] = 1
#             error[1].append('SiO2 must not be <40 or >49')
#     if TK-273.15 > 1500 or TK-273.15 < 600:
#         error[0] = 1
#         error[1].append('Temperature must be between 600 C and 1500 C.')
#==============================================================================
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Perform calculation if there is no error
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #initialize conditions
    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK, 10.0) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
    output = [out] #output records degassing path
    press = out[0]
    GH2O = out[5]
    GCO2 = out[6]
    routine = 1 #P is approximately known
    #print 'Degassing path progress: ' + str(round(100*1./steps,1)) + '%'
    if style == 0: #for open-system degassing
        WtH2Olast = WtH2O
        inc = WtH2O / 1000.
        Drop = PPMCO2 / steps

        for seq in range(1,steps,1): #Start for next loop to increment CO2 along degassing trend.
            loopcount = 0
            CarbonVFrac = 10.0 #Set to an arbitrarily high value to enter loop
        
            XX = 1.0
            PPMCO2o = PPMCO2 #Define previous CO2 concentration.
            PPMCO2 = PPMCO2 - Drop #Define new CO2 concentration.
            WTH2Oo = WtH2O #Define previous H2O concentration.
            
            while abs(CarbonVFrac - GCO2) >= 0.0003: #calculated by SatPress = that calculated by mass balance
                loopcount = loopcount + 1
                out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK, press) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                PPMCO2 = out[2]
                WtH2O = out[1]
                GCO2 = out[6]
                press = out[0]
                #routine[1] = press #reset pressure to the current output
                
                CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles of carbon in vapor by mass balance.
                WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles of water in vapor by mass balance.
                CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol fraction CO2 in vapor by mass balance.
                Xz = XX
                XX = copysignf(1.0, (-1 * (CarbonVFrac - GCO2)))
                if Xz + XX == 0:
                    inc = inc / 2.
                WtH2O = WtH2O + inc * XX #Use Newton's Method to cycle WtH2O until vapor composition
                if loopcount > 1000: #This if-then  re-initializes if there is trouble converging.
                    WtH2O = WtH2Olast
                    inc = WtH2O / 1000.
                    loopcount = 0
                  
            if WtH2Olast > WtH2O:
                inc = (WtH2Olast - WtH2O) / 30. #Define appropriate amount to increment H2O for next step.
            WtH2Olast = WtH2O
            output.append(out) #record of degassing path
            #print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'
    
    if style == 1: #for closed-system degassing
        
        WFCO2v = (GCO2 * 44.009) / (GH2O * 18.015 + GCO2 * 44.009) #Weight fraction CO2 in vapor.
        WtVap = excess / (1 - 0.01 * excess) #Grams excess vapor per 100 grams melt+ dissolved volatiles.
        
        MassVCO2 = WFCO2v * WtVap #Mass of CO2 in the excess vapor.
        MassVH2O = (1 - WFCO2v) * WtVap #Mass of H2O in the excess vapor.

        WTH2Oo = WtH2O + MassVH2O #Total grams H2O in system per 100 grams melt+dissolved volatiles at starting P.
        PPMCO2o = PPMCO2 + 10000 * MassVCO2 #Total grams CO2 in system per 100 grams melt+dissolved volatiles at starting P.
        WtH2Olast = WtH2O #Starting H2O dissolved in melt.

        Drop = PPMCO2 / steps #Divides CO2 into steps.
        inc = WtH2O / 1000. #Amount to vary WtH2O
        
        for seq in range(1,steps,1):
            loopcount = 0 #Keeps track of number of iterations. Re-initializes if loopcount >1000
            CarbonVFrac = 10.0 #Set to an arbitrarily high value to enter loop
            XX = 1.0
            PPMCO2 = PPMCO2 - Drop
            if PPMCO2 > 0:
                while abs(CarbonVFrac - GCO2) >= 0.0003:
                    loopcount = loopcount + 1
                    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK, press) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                    press = out[0]
                    PPMCO2 = out[2]
                    WtH2O = out[1]
                    GCO2 = out[6]
                    #routine[1] = press #reset pressure to the current output
                    
                    CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles CO2 in vapor by mass balance.
                    WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles H2O in vapor by mass balance.
                    CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol%CO2 in vapor
                    Xz = XX #Starts Newton's method to adjust WtH2O at given PPMCO2
                    XX = copysignf(1.0, (-1 * (CarbonVFrac - GCO2))) #until the vapor composition calculated by SatPress is equal
                    if Xz + XX == 0:
                        inc = inc / 2. #to that calculated by mass balance.
                    WtH2O = WtH2O + inc * XX
                    if loopcount > 1000: #This if-then  re-initializes if there is trouble converging
                        WtH2O = WtH2Olast
                        inc = WtH2O / 1000.
                        loopcount = 0
                        
                if WtH2Olast > WtH2O:
                    inc = (WtH2Olast - WtH2O) / 30.
                WtH2Olast = WtH2O
                
                output.append(out) #record of degassing path
                #print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'
    return output



