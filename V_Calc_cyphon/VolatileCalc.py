# -*- coding: utf-8 -*-

import numpy as np

def VolatileCalc(func,var1,var2):
    ###########################################################################
    #Acknowledgements
    ###########################################################################
    # This script is a direct port of the script developed in the following publication:
    #Newman, Sally, and Jacob B. Lowenstern. "VolatileCalc: a silicate melt–H2O–CO2 solution model written in Visual Basic for excel." Computers & Geosciences 28.5 (2002): 597-604.
    #Please cite the publication above for use of this script
    
    error = [0,[]] #tracks if an error occurs -[no(0)/yes(1),[list containing error messages]]
    
    
    ###########################################################################
    # Function parameter explanation
    ###########################################################################
    #func is a string that indicates the desired function
    if func.lower() in ['f','fug','fugacity']: #for calculations of vapor saturation pressure
        func = 0
    elif func.lower() in ['sp','saturation pressure']: #for calculations of vapor saturation pressure
        func = 1
    elif func.lower() in ['dp','degassing path']: #for calculations of degassing paths
        func = 2
    elif func.lower() in ['ib','isobar']: #for calculations of isobars
        func = 3
    elif func.lower() in ['ip','isopleth']: #for calculations of isopleths
        func = 4
    elif func.lower() in ['svp']: #for calculation of solubility vs pressure
        func = 5
    else:
        error[0] = 1
        error[1].append('Incorrect function selection. Must be sp (for saturation pressure), dp (for degassing path), ib (for isobar), or ip (for isopleth).')
        output = 'error'
    #In all cases except for calculation of fugacity:
        #var1 describes the composition of the system and must be one of 'basalt' or 'rhyolte'
    #Where fugacity is being calculated:
        #var1 is a float indicating T (in C)
    if func != 0: #comp =  0 (for basalt) or 1 (for rhyolite)
        if var1.lower() in ['basalt','bas','ba']: #for basaltic composition
            comp = 0
        elif var1.lower() in ['rhyolite','rhyol','rhy']: #for rhyolitic composition
            comp = 1
        else:
            error[0] = 1
            error[1].append('Incorrect composition selection. Must be basalt or rhyolite.')
    
    #if func = 'fug' (for fugacity calculation)
        #var1 = T (in C) - flaot
        #var2 = P (in MPa) - float

    #if func = 'sp' (for vapor sat pressure) and var1 = 'basalt', var2 must be list, which includes the following:
        #[H2O (in wt%),CO2 (in ppm),SiO2 (in wt%),Temperature (in C)]
        #SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'sp' (for vapor sat pressure) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[H2O (in wt%),CO2 (in ppm),Temperature (in C)]
        #Elements of the list must be float or integer
        
    #if func = 'dp' (for degassing path) and var1 = 'basalt', var2 must be list, which includes the following:
        #[H2Ostart (in wt%),CO2start (in ppm),SiO2 (in wt%),Temperature(in C), [style], steps]
        #SiO2 must not be less than 40 or greater than 49
        #the style list should be '[0]' for open system and '[1,#]' for closed system where # is the amount of excess vapor in weight percent (zero excess vapor should look like [1,0]), the first element should be an integer and the second is a float
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer
    #if func = 'dp' (for degassing path) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[H2Ostart (in wt%),CO2start (in ppm),Temperature(in C), [style], steps]
        #the style list should be '[0]' for open system and '[1,#]' for closed system where # is the amount of excess vapor in weight percent (zero excess vapor should look like [1,0]), the first element should be an integer and the second is a float
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer
    
    #if func = 'ib' (for isobar) and var1 = 'basalt', var2 must be list, which includes the following:
        #[P (in MPa), T (in C), SiO2 (in wt%), steps]
        #SiO2 must not be less than 40 or greater than 49
        #steps is an integer indicating the number of calculations
        #Other elements of the list must be float or integer
    #if func = 'ib' (for isobar) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[P (in MPa), T (in C), steps]
        #steps is an integer indicating the number of calculations
        #Elements of the list must be float or integer
        
    #if func = 'ip' (for isopleth) and var1 = 'basalt', var2 must be list, which includes the following:
        #[CO2 increment (in ppm), Molar percentage of H2O in fluid phase (in %), steps, SiO2 (in wt%), Temperature(in C)]
        #Notes: SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'ip' (for isopleth) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[CO2 increment (in ppm), Molar percentage of H2O in fluid phase (in %), steps, Temperature(in C)]
        #Elements of the list must be float or integer
        
    #if func = 'svp' (for solubility vs. pressure) and var1 = 'basalt', var2 must be list, which includes the following:
        #[variable, pressure interval (in MPa), SiO2 (in wt%), temperature (in C)]
        #variable is a string that must be 'h2o' or 'co2' (to indicate the volatile of interest)
        #Notes: SiO2 must not be less than 40 or greater than 49
        #Elements of the list must be float or integer
    #if func = 'svp' (for solubility vs. pressure) and var1 = 'rhyolite', var2 must be list, which includes the following:
        #[variable, pressure interval (in MPa), temperature (in C)]
        #variable is a string that must be 'h2o' or 'co2' (to indicate the volatile of interest)
        #Elements of the list must be float or integer
    ###########################################################################
    
    
    ###########################################################################
    # Internal functions
    ###########################################################################

    def FNA(TK):
        return (166800000 - 193080 * (TK - 273.15) + 186.4 * (TK - 273.15)**2 - 0.071288 * ((TK - 273.15)**3)) * 1.01325

    def FNB(TK):
        return 1.01325 * (73030000 - 71400 * (TK - 273.15) + 21.57 * (TK - 273.15)**2)

    def FNC(TK):
        R = 83.14321
        return 1.01325 * (np.exp(-11.071 + 5953 / TK - 2746000 / TK**2 + 464600000 / TK**3) * 0.5 * R * R * TK**2.5 / 1.02668 + 40123800)

    def FNF(V,TK,A,B,P):
        R = 83.14321
        return R * TK / (V - B) - A / ((V * V + B * V) * TK**0.5) - P

    def MRK(P,TK): #Redlich-Kwong routine to estimate endmember H2O and CO2 fugacities
        R = 83.14321
        B_1 = 14.6
        B_2 = 29.7
        for X_1 in range(2): #loops twice, once for each CO2 and H2O
            B = X_1 * B_1 + (1 - X_1) * B_2
            A = X_1**2 * FNA(TK) + 2 * X_1 * (1 - X_1) * FNC(TK) + (1 - X_1)**2 * FNB(TK)
            Temp2 = B + 5
            Q = 1
            Temp1 = 0
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
            G_1 = np.log(V / (V - B)) + B_1 / (V - B) - 2 * (X_1 * FNA(TK) + (1 - X_1) * FNC(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
            G_1 = G_1 + (np.log((V + B) / V) - B / (V + B)) * A * B_1 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
            G_1 = np.exp(G_1)
            G_2 = np.log(V / (V - B)) + B_2 / (V - B) - 2 * (X_1 * FNC(TK) + (1 - X_1) * FNB(TK)) * np.log((V + B) / V) / (R * TK**1.5 * B)
            G_2 = G_2 + (np.log((V + B) / V) - B / (V + B)) * A * B_2 / (R * TK**1.5 * B**2) - np.log(P * V / (R * TK))
            G_2 = np.exp(G_2)
            if X_1 == 0:
                fCO2o = G_2 * P #The fugacity of CO2
            if X_1 == 1:
                fH2Oo = G_1 * P #The fugacity of H2O
        return fCO2o, fH2Oo
    
    def CO2only(routine,comp,SiO2,PPMCO2,TK):
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #routine is a list, the first element is M
        #M = 0 (calc saturation pressure), 1 (equilibrium speciation), 2 (isobar calculation)
        #if M = 1 or 2, the second element in M is a starting pressure
        #comp = 0 (basalt), 1 (rhyolite)
        #SiO2 = SiO2 content in weight percent
        #PPMCO2 = CO2 content in ppm
        #TK = temperature in kelvin
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #local variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Determine if composition is basalt (0) or rhyolite (1)
        M = routine[0]
        if M != 0: #set pressure cases other than vapor sat pressure calculation
            press = routine[1]
        
        Z = 100. #This is the increment of CO2 (for case of M = 2, i.e. isobar calculation)
        WtH2O = 0 #Anhydrous run
        xH2Om = 0
        WtOHm = 0
        WtH2Om = 0
        
        H2Ocoeff = (-0.0000304 + 0.00000129 * SiO2) / 0.0000328
        CO2coeff = (0.0000087 - 0.0000001698 * SiO2) / 0.00000038
        
        
        if M != 2: #selects a pressure step and an approximate pressure to start calculations
            if PPMCO2 > 10 and PPMCO2 < 200:
                press = 200.
                Z = 1.
            elif PPMCO2 < 400:
                press = 500.
                Z = 100.
            elif PPMCO2 < 800:
                press = 750.
                Z = 100.
            elif PPMCO2 < 1600:
                press = 1500.
                Z = 100.
            elif PPMCO2 < 2400:
                press = 2500.
                Z = 100.
            else:
                press = 3000.
                Z = 100.
        elif M == 2: #M == 2 is for isobars, pressure is not variable
            if press <= 1000:
                PPMCO2 = 400.
            elif press <= 2000:
                PPMCO2 = 1200.
            elif press <= 3000:
                PPMCO2 = 2000.
            elif press <= 4000:
                PPMCO2 = 2800.
            elif press <= 5000:
                PPMCO2 = 3600.
            elif press <= 6000:
                PPMCO2 = 4400.
            else:
                PPMCO2 = 8000.
        
        #initialize Y for Newton's method
        Y = 1
        GH2O = 2
        GCO2 = 0
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #isobar loop
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        while abs(GH2O + GCO2 - 1) >= 0.0001:
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate xCO2m
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if GH2O == 2 or M==2: #Run this part of the routine always in the first loop, then repeat for isobar calculation
                SNCO2 = PPMCO2 * 0.0001 / 44.009 #Relative moles of CO2
                SNH2O = WtH2O / 18.015 #Relative moles of H2O: zero for anhydrous
                if comp == 0: #if basaltic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                else: #if rhyolitic
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 #relative moles of oxygen
                XH2O = SNH2O / (SNH2O + SNO) #Mole Fraction H2O = zero for anhydrous
                XO = 1 - XH2O #Mole fraction oxygen = 1 for anhydrous
                xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O)
                if comp == 0: #if basaltic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
                else: #if rhyolitic
                    xb = SNH2O / (SNH2O + (100 - WtH2O - PPMCO2 * 0.0001) / 32.5 + (PPMCO2 * 0.0001 / 44.009)) #xb is SNH2O over (SNH2O+SNO+SNCO2)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #Calculate gas
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            fCO2o,fH2Oo = MRK(press,TK)
            #
            #GH2O and GCO2 are the mol fractions of H2O and CO2 in the gas phase.  They are dependent on XCO2m and XH2Om
            #as well as T and P and endmember H2O and CO2 fugacities.
            #Standard state for basalt is 1 bar and 1200C.  All values, plus Delta V and Delta H from Dixon et al. (1995).
            #Solubilities at standard state are 0.5 ppm CO2 and 0.11 wt.% H2O for a basalt with 49% SiO2. These values
            #are equivalent to mol fractions of .00000038 and .0000328, respectively.
            #Solubilities vary accoding to wt.% SiO2 in accord with H2Ocoeff and CO2coeff derived by Dixon (1997) and
            #slightly modified so that a basalt with 49% SiO2 is equivalent to the tholeiite in the original (1995) model.
            #
            #If GH2O+GCO2 sum to > 1, then pressure must be decreased or either Wt.H2O or PPMCO2 need to be decreased.
            #
            
            if comp == 0: #if composition is basaltic
                GH2O = (xH2Om / 0.0000328 / H2Ocoeff) * (1 / fH2Oo) / np.exp(-12 * (press - 1) / (41.84 * 1.9872 * TK))
                GCO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(-23 * (press - 1) / (41.84 * 1.9872 * TK))
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
                GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(-5 * (press - 799) / (41.84 * 1.9872 * TK) + 4420 / 1.9872 * (1 / TK - 1 / 1123.16))
                GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(-28 * (press - 750) / (41.84 * 1.9872 * TK) + 4861 * (1 / TK - 1 / 1123.16) / 1.9872)
            
            #redifine variables
            if abs(GH2O + GCO2 - 1) >= 0.0001:
                Yo = Y
                Y = np.sign(-1. * (GH2O + GCO2 - 1))
                if Y + Yo == 0:
                    Z = Z/2.
                if M == 2:
                    PPMCO2 = PPMCO2 + Z * Y / 2.
                else:
                    press = press - Z * Y / 2.
        
        return [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
            
        
    def SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK):
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
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #local variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #If vapor saturation pressure calculation, routine[0] = 1
        M = routine[0]
        if M != 0: #set pressure cases other than vapor sat pressure calculation
            press = routine[1]
        
        Z = 10. #This is the increment of pressure, default is 10 bar increment
        
        H2Ocoeff = (-0.0000304 + 0.00000129 * SiO2) / 0.0000328
        CO2coeff = (0.0000087 - 0.0000001698 * SiO2) / 0.00000038
        
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
        Y = 1
        GH2O = 2
        GCO2 = 0
        
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
                xb = 0
                XOHm = 0
                if comp == 0: #if the composition is basaltic
                    xH2Om = np.exp(-5.827) * WtH2O**1.855
                    SNCO2 = PPMCO2 * 0.0001 / 44.009 #relative moles of CO2
                    SNH2O = WtH2O / 18.015 #relative moles of H2O
                    SNO = (100 - WtH2O - PPMCO2 * 0.0001) / 36.594 #relative moles of oxygen
                    xCO2m = SNCO2 / (SNCO2 + SNO + SNH2O) #mol fraction of CO2 relative to H2O+CO2+O
                    
                else: #if the composition is rhyolitic
                    xH2Om = np.exp(-5.526) * WtH2O**1.977
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
                derath = 1
                if comp == 0: #if composition is basaltic
                    while abs(derath) >= 0.00001: #loop to define XOHm
                        fx = 9.143 - 3.295 * (XOHm - 1) - 2 * 6.019 * (XO - XOHm) - 2 * 0.572 * (XH2O - XOHm) + np.log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
                        fxp = -3.295 + 2 * (6.019 + 0.572) + (2 * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm) + 0.5 * XOHm * ((XH2O - 0.5 * XOHm) + (XO - 0.5 * XOHm))) / (XOHm * (XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm))
                        derath = fx / fxp
                        XOHm = XOHm - derath
                else: #if composition is rhyolitic
                    while abs(derath)>= 0.00001: #loop to define XOHm
                        fx = 9.345 - 4.304 * (XOHm - 1) - 2 * 6.277 * (XO - XOHm) - 2 * 2.328 * (XH2O - XOHm) + np.log(XOHm**2 / ((XH2O - 0.5 * XOHm) * (XO - 0.5 * XOHm)))
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
                GH2O = (xH2Om / 0.0000328 / H2Ocoeff) * (1 / fH2Oo) / np.exp(-12 * (press - 1) / (41.84 * 1.9872 * TK))
                GCO2 = (xCO2m / 0.00000038 / CO2coeff) * (1 / fCO2o) / np.exp(-23 * (press - 1) / (41.84 * 1.9872 * TK))
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
                GH2O = (xH2Om / 0.0323) * (712.2 / fH2Oo) / np.exp(-5 * (press - 799) / (41.84 * 1.9872 * TK) + 4420 / 1.9872 * (1 / TK - 1 / 1123.16))
                GCO2 = (xCO2m / 0.000397) * (901.6 / fCO2o) / np.exp(-28 * (press - 750) / (41.84 * 1.9872 * TK) + 4861 * (1 / TK - 1 / 1123.16) / 1.9872)
            
            #Redefine variables
            if abs(GH2O + GCO2 - 1) >= 0.0001:
                Yo = Y
                Y = np.sign(-1. * (GH2O + GCO2 - 1))
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
        
        return [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]

    ###########################################################################
    # Fugacity
    ###########################################################################
    if error[0] == 0 and func == 0:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TK = var1+273.15
        P = float(var2)*10
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation and test for error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        if P >= 200 and TK >= 450+273.15:
            fCO2o,fH2Oo = MRK(P,TK)
            output = [fH2Oo, fCO2o]
        else:
            error = [1,'This modified Redlich-Kwong is not recommended for pressures < 200 bars and temperatures <450°C']

    ###########################################################################
    # Saturation pressure
    ###########################################################################
    if error[0] == 0 and func == 1:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WtH2O = float(var2[0])
        PPMCO2 = float(var2[1])
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
        else:
            SiO2 = 0
            TK = var2[2]+273.15

        routine = [0] #routine = [0] means that SatPress function varies pressure
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if WtH2O == 0 and PPMCO2 <= 10:
            error[0] = 1
            error[1].append('For anhydrous runs, VolatileCalc needs a CO2  concentration > 10 ppm')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            if WtH2O == 0:
                output = CO2only(routine,comp,SiO2,PPMCO2,TK) #output = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
            else:
                output = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #output = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
        
        
    ###########################################################################
    # Degassing path
    ###########################################################################
    elif error[0] == 0 and func == 2:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        WtH2O = float(var2[0])
        PPMCO2 = float(var2[1])
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
            style = int(var2[4][0]) #var[3] is a list, the first element is 0 = open, 1 = closed
            if len(var2[4]) > 1: #if closed there is a second element that contains wt% of excess vapor
                excess = float(var2[4][1])
            else: #Default excess vapor is 0 wt% (only applicable when considering closed system degassing)
                excess = 0
            steps = int(var2[5])
        else:
            SiO2 = 0
            TK = var2[2]+273.15
            style = int(var2[3][0]) #var[3] is a list, the first element is 0 = open, 1 = closed
            if len(var2[3]) > 1: #if closed there is a second element that contains wt% of excess vapor
                excess = float(var2[3][1])
            else: #Default excess vapor is 0 wt% (only applicable when considering closed system degassing)
                excess = 0
            steps = int(var2[4])

        routine = [0] #routine = [0] means that SatPress function varies pressure

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            #initialize conditions
            out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
            output = [out] #output records degassing path
            press = out[0]
            GH2O = out[5]
            GCO2 = out[6]
            routine = [1,press] #P is approximately known
            print 'Degassing path progress: ' + str(round(100*1./steps,1)) + '%'
            if style == 0: #for open-system degassing
                WtH2Olast = WtH2O
                inc = WtH2O / 1000.
                Drop = PPMCO2 / steps
                for seq in range(1,steps,1): #Start for next loop to increment CO2 along degassing trend.
                    loopcount = 0
                    CarbonVFrac = 10 #Set to an arbitrarily high value to enter loop
                
                    XX = 1
                    PPMCO2o = PPMCO2 #Define previous CO2 concentration.
                    PPMCO2 = PPMCO2 - Drop #Define new CO2 concentration.
                    WTH2Oo = WtH2O #Define previous H2O concentration.
                    
                    while abs(CarbonVFrac - GCO2) >= 0.0003: #calculated by SatPress = that calculated by mass balance
                        loopcount = loopcount + 1
                        out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                        PPMCO2 = out[2]
                        WtH2O = out[1]
                        GCO2 = out[6]
                        press = out[0]
                        routine[1] = press #reset pressure to the current output
                        
                        CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles of carbon in vapor by mass balance.
                        WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles of water in vapor by mass balance.
                        CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol fraction CO2 in vapor by mass balance.
                        Xz = XX
                        XX = np.sign(-1 * (CarbonVFrac - GCO2))
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
                    print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'
            
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
                    CarbonVFrac = 10 #Set to an arbitrarily high value to enter loop
                    XX = 1
                    PPMCO2 = PPMCO2 - Drop
                    if PPMCO2 > 0:
                        while abs(CarbonVFrac - GCO2) >= 0.0003:
                            loopcount = loopcount + 1
                            out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #out = [press, WtH2O, PPMCO2, WtH2Om, WtOHm, GH2O, GCO2]
                            PPMCO2 = out[2]
                            WtH2O = out[1]
                            GCO2 = out[6]
                            routine[1] = press #reset pressure to the current output
                            
                            CarbonVapor = (PPMCO2o - PPMCO2) / 44.009 #Calculate moles CO2 in vapor by mass balance.
                            WaterVapor = (WTH2Oo - WtH2O) * 10000 / 18.02 #Calculate moles H2O in vapor by mass balance.
                            CarbonVFrac = CarbonVapor / (CarbonVapor + WaterVapor) #Calculate mol%CO2 in vapor
                            Xz = XX #Starts Newton's method to adjust WtH2O at given PPMCO2
                            XX = np.sign(-1 * (CarbonVFrac - GCO2)) #until the vapor composition calculated by SatPress is equal
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
                        print 'Degassing path progress: ' + str(round(100*(seq+1.)/steps,1)) + '%'
    
    
    ###########################################################################
    # Isobar calculation
    ###########################################################################
    elif error[0] == 0 and func == 3:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        P = var2[0]*10.
        TK = var2[1]+273.15
        if comp == 0:
            SiO2 = float(var2[2])
            steps = int(var2[3])
        else:
            SiO2 = 0
            steps = int(var2[2])
        
        routine = [2,P] #Signals Saturation Pressure function that P is known
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if P > 5000:
            error[1].append('VolatileCalc is not recommended for P > 500 MPa.')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = []
            final = CO2only(routine,comp,SiO2,0,TK)
            MaxCO2 = final[2]
            if P < 500:
                WtH2O = 2.
            elif P < 1000:
                WtH2O = 3.
            elif P < 2000:
                WtH2O = 5.
            elif P < 3000:
                WtH2O = 7.
            elif P < 4000:
                WtH2O = 8.
            else:
                WtH2O = 9.
            
            addCO2 = MaxCO2/(steps-1)
            PPMCO2 = 0
            
            for W in range(steps):
                if W != steps-1:
                    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                    WtH2O = out[1]
                    PPMCO2 += addCO2
                else:
                    out = final
                output.append(out)
    
    ###########################################################################
    # Isopleth calculation
    ###########################################################################
    elif error[0] == 0 and func == 4:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CO2inc = int(var2[0])
        GH = float(var2[1]) / 100. #Molar percent of H2O in fluid converted to fraction
        steps = int(var2[2])
        if comp == 0:
            SiO2 = float(var2[3])
            TK = var2[4]+273.15
        else:
            SiO2 = 0
            TK = var2[3]+273.15

        routine = [0] #routine = [0] means that SatPress function varies pressure
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = [[0, 0, 0, 0, 0, GH, 1 - GH]] #List to store output, first point is at the origin
            WtH2O = 1.0 #Guess 1wt% H2O to start
            PPMCO2 = CO2inc #First calc is 1 CO2 step
            for i in range(steps-1):
                YY = 0 #YY tells whether H2O overestimated or underestimated.
                ZZ = 0.2 #Initial amount to change WtH2O.
                out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                press = out[0]
                GH2O = out[-2]
                routine = [1,press]
                
                #The following loop adjusts WtH2O given a constant PPMCO2, T and approximate P.
                while abs(GH2O - GH) >= 0.0001: #Stop when calculated vapor comp is same as original requested vapor comp.
                    Yi = YY
                    YY = np.sign(-1 * (GH2O - GH)) #Is calc. vapor comp > or < desired vapor comp?
                    if YY + Yi == 0: #Newton's method.
                        ZZ = ZZ / 2.
                    WtH2O = WtH2O + ZZ * YY / 2. #Adjust water in melt.
                    out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK) #Calculate new vapor comp.
                    press = out[0]
                    GH2O = out[-2]
                    routine = [1,press]
                    
                output.append(out)
                
                PPMCO2 += CO2inc #Prepare next CO2 concentration

    ###########################################################################
    # Solubility vs pressure calculation
    ###########################################################################
    elif error[0] == 0 and func == 5:
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Unpack input variables
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        var = var2[0].lower()
        P_interval = int(var2[1]*10)
        if comp == 0:
            SiO2 = float(var2[2])
            TK = var2[3]+273.15
        else:
            SiO2 = 0
            TK = var2[2]+273.15

        routine = [2,5000] #routine = [2] means that in the sat press functions, P is known
    
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Test input variables for errors
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if comp == 0:
            if SiO2 < 40 or SiO2 > 49:
                error[0] = 1
                error[1].append('SiO2 must not be <40 or >49')
        if TK-273.15 > 1500 or TK-273.15 < 600:
            error[0] = 1
            error[1].append('Temperature must be between 600 C and 1500 C.')
        if var == 'h2o':
            GH = 1
            WtH2O = 1 #initial guess at H2O
            PPMCO2 = 0
        elif var == 'co2':
            GH = 0
            PPMCO2 = 100 #initial guess at CO2
        else:
            error[0] = 1
            error[1].append('For solubility vs. pressure calculations, you must request H2O or CO2.')
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Perform calculation if there is no error
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if error[0] == 0:
            output = [[0, 0, 0, 0, 0, GH, 1 - GH]] #List to store output, first point is at the origin
            steps = 5000/P_interval #Set steps
            P_start = P_interval
            P_stop = P_interval * steps
            for P in range(P_start,P_stop+P_interval,P_interval):
                routine[1] = float(P)
                if GH == 1:
                     out = SatPress(routine,comp,SiO2,WtH2O,PPMCO2,TK)
                     out[-2] = 1
                else:
                    out = CO2only(routine,comp,SiO2,PPMCO2,TK)
                    out[-1] = 1
                output.append(out)

    ###########################################################################
    # Prepare output
    ###########################################################################
    
    if error[0] == 1:
        output = 'Error'
    elif func == 0: #Fugacity
        output[0] /= 10 #Convert from bars to MPa
        output[1] /= 10 #Convert from bars to MPa
    elif func == 1: #Vapor saturation pressure
        output[0] /= 10 #Convert from bars to MPa
        output[-2] *= 100 #Convert H2Ov from fraction to percent
        output[-1] *= 100 #Convert CO2v from fraction to percent
        if output[0] > 500:
            error[1].append('VolatileCalc is not recommended for P > 500 MPa.')
    elif func in [2,3,4,5]: #Degassing path, isobar
        for i in range(len(output)):
            output[i][0] /= 10 #Convert from bars to MPa
            output[i][-2] *= 100 #Convert H2Ov from fraction to percent
            output[i][-1] *= 100 #Convert CO2v from fraction to percent
            if output[i][0] > 500 and 'VolatileCalc is not recommended for P > 500 MPa.' not in error[1]:
                error[1].append('VolatileCalc is not recommended for P > 500 MPa.')
    
    return output, error[1]