# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 14:39:42 2025

@author: Adhikari
"""

import numpy as np
import os
from multiprocessing import Pool, Process
from subprocess import Popen, PIPE, CalledProcessError, run
import math
import subprocess
import pandas as pd
import shutil 
import matplotlib.pyplot as plt

def param_blocks(param_values, mc_cores = None):
    if mc_cores == None:
        mc_cores = os.cpu_count()
    if mc_cores > os.cpu_count():
        mc_cores = os.cpu_count()
    
    param_blocks = []
    work_div = param_values.shape[0]/mc_cores
    if work_div - int(work_div) >= 0.5:
        mc_index = int(math.ceil(param_values.shape[0]/mc_cores))
    elif work_div - int(work_div) < 0.5:
        mc_index = int(math.floor(param_values.shape[0]/mc_cores))
    start = 0
    for i in range(mc_cores):
        if i < mc_cores-1:
            end = start + mc_index
            #print(start,end)
            #param_small_index.append([start, end])
            param_blocks.append(param_values[start:end])
            start = end
        else:                                                                       #bypass if params cannot be divided into 16 equal parts
            end = param_values.shape[0]
            #print(start,end)
            #param_small_index.append([start:,:])
            param_blocks.append(param_values[start:])
    return param_blocks

def irrigation_threshold_calculator(rootdepth, nfk_drop):
    
    #MvG parameters, and soil pedo transfer function parameters
    h_pwp = -10**4.2
    h_fc = -10**1.8
    ores = np.array([0.0387, 0.0375, 0.0396, 0.04, 0.0473])
    osat = np.array([0.4376, 0.3446, 0.3576, 0.3603, 0.3424])
    alpha = np.array([0.0459, 0.0463, 0.0423, 0.0339, 0.0348])
    n = np.array([1.7434, 1.9512, 1.566, 1.4164, 3.3022])
    m = np.array([1 - 1 / val for val in n])
    pwp = ores + (osat - ores) / ((1 + np.abs(alpha * h_pwp)**n)**m)
    fc = ores + (osat - ores) / ((1 + np.abs(alpha * h_fc)**n)**m)
    
    
    #soil horizon depths
    horizons = [
        [0, 30],
        [30, 50],
        [50, 60],
        [60, 95],
        [95, 200]]
    #weight = depth of mature root within the soil horizon 
    weight_list = [0] * len(ores)
    for i in range(len(horizons)):
        low = horizons[i][0]
        high = horizons[i][1]
        if rootdepth >= high:
            weight_list[i] = high - low
        elif rootdepth > low and rootdepth < high:
            weight_list[i] = rootdepth - low
        elif rootdepth < low:       #else statement also works here, but more logical control and clarity 
            weight_list[i] = 0
    weight = np.array(weight_list)
    # --- Threshold moisture calculation ---
    oall = pwp + nfk_drop/100 * (fc-pwp)
    wc_thres = (oall * weight).sum()*10 #cm to mm
    return wc_thres

def write_management_nsga_II(dir_, parameter, calibration):
    initialize = f'''(description "Simulation for use in Hamerstorf.")
   (path &old
 "C:/Program Files/Daisy 6.45/sample" 
 "C:/Program Files/Daisy 6.45/lib")
   ;; Using standard parameterizations 
   (input file "tillage.dai")
   (input file "crop.dai") 
   (input file "log.dai")
   (input file "fertilizer.dai")'''
   
    soil_horizons = f''';; creating soil horizons
   (defhorizon Ap FAO3
 "Hamerstorf top soil."
 (clay 2.0 [%])
 (silt 16.0 [%])
 (sand 82.0 [%])
 (humus 1.12 [%])
 (hydraulic M_vG
   (K_sat 16.33 [cm/h])
   (Theta_res 3.87 [%])
   (Theta_sat 43.76 [%])
   (alpha 0.0454 [cm^-1])
   (n 1.7434))
 (dry_bulk_density 1.394 [g/cm^3]))
 
   (defhorizon B Ap
   "Hamerstorf Second Layer"
   (clay 2.0 [%])
   (silt 16.0 [%])
   (sand 82.0 [%])
   (humus 0.29 [%])
   (hydraulic M_vG
   (K_sat 4.75 [cm/h])
   (Theta_res 3.75 [%])
   (Theta_sat 34.46 [%])
   (alpha 0.0463 [cm^-1])
   (n 1.9512))
   (dry_bulk_density 1.624 [g/cm^3]))
   
   (defhorizon C Ap
   "Hamerstorf Third Layer"
   (clay 7.0 [%])
   (silt 18.0 [%])
   (sand 75.0 [%])
   (humus 0.01 [%])
   (hydraulic M_vG
   (K_sat 2.3 [cm/h])
   (Theta_res 3.96 [%])
   (Theta_sat 35.76 [%])
   (alpha 0.0423 [cm^-1])
   (n 1.566))
   (dry_bulk_density 1.597 [g/cm^3]))

   (defhorizon D Ap
   "Hamerstorf Fourth Layer"
   (clay 10.0 [%])
   (silt 25.0 [%])
   (sand 65.0 [%])
   (humus 0.01 [%])
   (hydraulic M_vG
   (K_sat 1.3 [cm/h])
   (Theta_res 4.0 [%])
   (Theta_sat 36.03 [%])
   (alpha 0.0339 [cm^-1])
   (n 1.4164))
   (dry_bulk_density 1.575 [g/cm^3]))

   (defhorizon E Ap
   "Hamerstorf Fifth Layer"
   (clay 1.0 [%])
   (silt 5.0 [%])
   (sand 94.0 [%])
   (humus 0.01 [%])
   (hydraulic M_vG
   (K_sat 22.21 [cm/h])
   (Theta_res 4.73 [%])
   (Theta_sat 34.24 [%])
   (alpha 0.0348 [cm^-1])
   (n 3.3022))
   (dry_bulk_density 1.64 [g/cm^3]))'''
   
   
    soil_column = f'''(defcolumn Hamerstorf default
 "Data collected by Hamerstorf."
 (Surface 
  (EpFactor {0.6} [])
  (albedo_dry {0.15} [])
  (albedo_wet {0.08} []))
 
 (Bioclimate default (pet PM))
 (Soil (horizons (-30 [cm] Ap) (-50 [cm] B) (-60 [cm] C) (-95 [cm] D) (-200 [cm] E))
       (MaxRootingDepth 130 [cm]))
 (Movement vertical(Geometry (zplus -2.5 -5 -10 -15 -20 -30 -40 -50 -60 -70 -80 -90  -100 -125 -150 -180 -200 [cm])))
 (Groundwater deep))'''
   

    crop_normal =f''';Creating crop parameters  
 (defcrop "WinterWheat_normal" "Winter Wheat"
   
  (Devel default 
(EmrTSum {parameter[0]}) (DSRate1 {parameter[1]}) (DSRate2 {parameter[2]})
(DS_Emr 0.01)
(TempEff1 (-10.0  0.01) (0.00  0.01) (20.0  0.90) (25.0  1.00) (35.0  1.20))
(PhotEff1 (10.0  0.29) (11.0  0.55) (12.0  0.75) (13.0  0.89) (14.0  1.00) (15.0  1.08) (16.0  1.14) (17.0  1.18) (24.0  1.18))
(TempEff2 ( 0.0  0.00) (10.0  0.14) (15.0  0.66) (20.0  1.00) (25.0  1.23))
   )
   
  (Vernal default 
(DSLim 0.25) (TaLim 5) (TaSum -50)
   )

  (LeafPhot original 
(Fm {parameter[3]}) 
(Qeff   {parameter[4]})
(TempEff (-20.0  0.00) (4.00  0.00) (10.0  1.00) (25.0  1.00) (35.0  0.01) (50.0  0.00))
   )

  (Seed release
(DM_fraction {87} [%])
(C_fraction {45} [%])
(N_fraction {2} [%])
(rate {0.4} [d^-1])
   )
   
  (Canopy
(SpLAI    {parameter[5]})
(PARref   0.06)
(PARext   0.6)
(PARrel   0.05)
(EPext    0.5)
(LAIDist0 0.00 0.00 1.00)
(LAIDist1 0.00 0.2 0.9)
(HvsDS   (0.00    1)    (1.00  80)    (2.00 100))
   )
  
   (Root (DptEmr     10)
(PenPar1    {parameter[6]})
(PenPar2     {parameter[7]})
(MaxPen      130)
(Rad        0.005)
(h_wp       -15000)
(Rxylem     10)
(MxNH4Up     2.5E-0007)
(MxNO3Up     2.5E-0007)
   )
   
   (Partit 
 (Root (0.0 {parameter[8]}) (0.3 {parameter[9]}) (0.5 0.4) (0.8 0.16) (0.95 0.04) (1.0 0.00) (2.0 0.00))
 (Leaf (0.0 {[parameter[10] if parameter[10] + 0.28 < 1 else 1-0.28 for i in range(1)][0]}) (0.3 {[parameter[11] if parameter[11] + 0.42 < 1 else 1-0.42 for i in range(1)][0]}) (0.5 {[parameter[12] if parameter[12] + 0.47 < 1 else 1-0.47 for i in range(1)][0]}) (0.8 0.3) (0.95 0.25) (1.1 0.00) (2.0 0.00))
 (Stem (0.0 0.28) (0.3 0.42) (0.5 0.47) (0.8 0.25) (0.95 0.62) (1.0 0.54) (1.1 0.38) (1.5 0.00) (2.0 0.00))         
 (NNI_crit 1.40 []); RATJEN & KAGE: 1.38
 (NNI_inc 0.65 []); RATJEN & KAGE: 0.60
 (RSR (0.0 1.1) (0.05 1) (0.3 0.5) (1.0 0.5) (2.0 0.25))
   )
   
   
   (Prod (NCrop 0.4)
(E_Root     0.69)
(E_Leaf     {parameter[13]})
(E_Stem     0.66)
(E_SOrg     0.7)
(r_Root     0.015)
(r_Leaf     0.016)
(r_Stem     0.01)
(r_SOrg     {parameter[14]})
(ShldResC   0.4)
(ReMobilDS  1.1)
(ReMobilRt  0.03)
(Large_RtDR 0.05)
(ExfoliationFac 0.7)
(LfDR (0.0 0.00) (0.3 0.00) (0.5 0.015) (1.0 {parameter[15]}) (1.5 {parameter[16]}) (2.0 0.08))
(RtDR (0.0 0.00) (0.5 0.004) (1.0 0.02) (2.0 0.03))
)
   
   (CrpN 
   (PtRootCnc (0.2 0.025) (0.5 0.026) (0.8 0.018) (1.20 0.015) (2.00 0.014))
   (CrRootCnc (0.2 0.017) (0.5 0.019) (0.8 0.014) (1.20 0.012) (2.00 0.009))
   (NfRootCnc (0.2 0.009) (1.20 0.003) (2.00 0.003))    
   (PtLeafCnc (0.2 0.051) (0.5 0.048) (0.8 0.044) (1.00 0.038) (2.00 0.021))
   (CrLeafCnc (0.2 0.039) (0.5 0.036) (1.00 0.025) (2.00 0.014))
   (NfLeafCnc (0.2 0.01) (0.5 0.01) (1.00 0.006) (2.00 0.005))
   (PtStemCnc (0.2 0.038) (0.5 0.035) (0.8 0.02) (1.00 0.015) (2.00 0.011))
   (CrStemCnc (0.2 0.023) (0.5 0.013) (0.8 0.008) (1.00 0.007) (2.00 0.006))
   (NfStemCnc (0.2 0.01) (1.00 0.002) (2.00 0.002))
   (PtSOrgCnc (1.00 0.017) (2.00 0.0175))
   (CrSOrgCnc (1.00 0.0125) (2.00 0.0125))
   (NfSOrgCnc (1.00 0.0005) (2.00 0.0005))
   (TLLeafEff (0 0.9) (2 0.9))
   (TLRootEff (0 0.1) (2 0.1))
   )
 
(Harvest (EconomicYield_W 0.8) (EconomicYield_N 0.94))
 )'''

    crop_wet = f''';Creating crop parameters  
 (defcrop "WinterWheat_wet" "Winter Wheat"
    (Devel default 
(EmrTSum {parameter[18]}) (DSRate1 {parameter[19]}) (DSRate2 {parameter[20]})
(DS_Emr 0.01)
(TempEff1 (-10.0  0.01) (0.00  0.01) (20.0  0.90) (25.0  1.00) (35.0  1.20))
(PhotEff1 (10.0  0.29) (11.0  0.55) (12.0  0.75) (13.0  0.89) (14.0  1.00) (15.0  1.08) (16.0  1.14) (17.0  1.18) (24.0  1.18))
(TempEff2 ( 0.0  0.00) (10.0  0.14) (15.0  0.66) (20.0  1.00) (25.0  1.23))
   )
   
  (Vernal default 
(DSLim 0.25) (TaLim 5) (TaSum -50)
   )

  (LeafPhot original 
(Fm {parameter[21]}) 
(Qeff   {parameter[22]})
(TempEff (-20.0  0.00) (4.00  0.00) (10.0  1.00) (25.0  1.00) (35.0  0.01) (50.0  0.00))
   )

  (Seed release
(DM_fraction {87} [%])
(C_fraction {45} [%])
(N_fraction {2} [%])
(rate {0.4} [d^-1])
   )
   
  (Canopy
(SpLAI    {parameter[23]})
(PARref   0.06)
(PARext   0.6)
(PARrel   0.05)
(EPext    0.5)
(LAIDist0 0.00 0.00 1.00)
(LAIDist1 0.00 0.2 0.9)
(HvsDS   (0.00    1)    (1.00  80)    (2.00 100))
   )
  
   (Root (DptEmr     10)
(PenPar1    {parameter[24]})
(PenPar2     {parameter[25]})
(MaxPen      130)
(Rad        0.005)
(h_wp       -15000)
(Rxylem     10)
(MxNH4Up     2.5E-0007)
(MxNO3Up     2.5E-0007)
   )
   
   (Partit 
 (Root (0.0 {parameter[26]}) (0.3 {parameter[27]}) (0.5 0.4) (0.8 0.16) (0.95 0.04) (1.0 0.00) (2.0 0.00))
 (Leaf (0.0 {[parameter[27] if parameter[27] + 0.28 < 1 else 1-0.28 for i in range(1)][0]}) (0.3 {[parameter[28] if parameter[28] + 0.42 < 1 else 1-0.42 for i in range(1)][0]}) (0.5 {[parameter[29] if parameter[29] + 0.47 < 1 else 1-0.47 for i in range(1)][0]}) (0.8 0.3) (0.95 0.25) (1.1 0.00) (2.0 0.00))
 (Stem (0.0 0.28) (0.3 0.42) (0.5 0.47) (0.8 0.25) (0.95 0.62) (1.0 0.54) (1.1 0.38) (1.5 0.00) (2.0 0.00))         
 (NNI_crit 1.40 []); RATJEN & KAGE: 1.38
 (NNI_inc 0.65 []); RATJEN & KAGE: 0.60
 (RSR (0.0 1.1) (0.05 1) (0.3 0.5) (1.0 0.5) (2.0 0.25))
   )
   
   
   (Prod (NCrop 0.4)
(E_Root     0.69)
(E_Leaf     {parameter[31]})
(E_Stem     0.66)
(E_SOrg     0.7)
(r_Root     0.015)
(r_Leaf     0.016)
(r_Stem     0.01)
(r_SOrg     {parameter[32]})
(ShldResC   0.4)
(ReMobilDS  1.1)
(ReMobilRt  0.03)
(Large_RtDR 0.05)
(ExfoliationFac 0.7)
(LfDR (0.0 0.00) (0.3 0.00) (0.5 0.015) (1.0 {parameter[33]}) (1.5 {parameter[34]}) (2.0 0.08))
(RtDR (0.0 0.00) (0.5 0.004) (1.0 0.02) (2.0 0.03))
)
   
   (CrpN 
   (PtRootCnc (0.2 0.025) (0.5 0.026) (0.8 0.018) (1.20 0.015) (2.00 0.014))
   (CrRootCnc (0.2 0.017) (0.5 0.019) (0.8 0.014) (1.20 0.012) (2.00 0.009))
   (NfRootCnc (0.2 0.009) (1.20 0.003) (2.00 0.003))    
   (PtLeafCnc (0.2 0.051) (0.5 0.048) (0.8 0.044) (1.00 0.038) (2.00 0.021))
   (CrLeafCnc (0.2 0.039) (0.5 0.036) (1.00 0.025) (2.00 0.014))
   (NfLeafCnc (0.2 0.01) (0.5 0.01) (1.00 0.006) (2.00 0.005))
   (PtStemCnc (0.2 0.038) (0.5 0.035) (0.8 0.02) (1.00 0.015) (2.00 0.011))
   (CrStemCnc (0.2 0.023) (0.5 0.013) (0.8 0.008) (1.00 0.007) (2.00 0.006))
   (NfStemCnc (0.2 0.01) (1.00 0.002) (2.00 0.002))
   (PtSOrgCnc (1.00 0.017) (2.00 0.0175))
   (CrSOrgCnc (1.00 0.0125) (2.00 0.0125))
   (NfSOrgCnc (1.00 0.0005) (2.00 0.0005))
   (TLLeafEff (0 0.9) (2 0.9))
   (TLRootEff (0 0.1) (2 0.1))
   )
 
(Harvest (EconomicYield_W 0.8) (EconomicYield_N 0.94))
)'''

    crop_dry = f''';Creating crop parameters  
 (defcrop "WinterWheat_dry" "Winter Wheat"
   
  (Devel default 
(EmrTSum {parameter[36]}) (DSRate1 {parameter[37]}) (DSRate2 {parameter[38]})
(DS_Emr 0.01)
(TempEff1 (-10.0  0.01) (0.00  0.01) (20.0  0.90) (25.0  1.00) (35.0  1.20))
(PhotEff1 (10.0  0.29) (11.0  0.55) (12.0  0.75) (13.0  0.89) (14.0  1.00) (15.0  1.08) (16.0  1.14) (17.0  1.18) (24.0  1.18))
(TempEff2 ( 0.0  0.00) (10.0  0.14) (15.0  0.66) (20.0  1.00) (25.0  1.23))
   )
   
  (Vernal default 
(DSLim 0.25) (TaLim 5) (TaSum -50)
   )

  (LeafPhot original 
(Fm {parameter[39]}) 
(Qeff   {parameter[40]})
(TempEff (-20.0  0.00) (4.00  0.00) (10.0  1.00) (25.0  1.00) (35.0  0.01) (50.0  0.00))
   )

  (Seed release
(DM_fraction {87} [%])
(C_fraction {45} [%])
(N_fraction {2} [%])
(rate {0.4} [d^-1])
   )
   
  (Canopy
(SpLAI    {parameter[41]})
(PARref   0.06)
(PARext   0.6)
(PARrel   0.05)
(EPext    0.5)
(LAIDist0 0.00 0.00 1.00)
(LAIDist1 0.00 0.2 0.9)
(HvsDS   (0.00    1)    (1.00  80)    (2.00 100))
   )
  
   (Root (DptEmr     10)
(PenPar1    {parameter[42]})
(PenPar2     {parameter[43]})
(MaxPen      130)
(Rad        0.005)
(h_wp       -15000)
(Rxylem     10)
(MxNH4Up     2.5E-0007)
(MxNO3Up     2.5E-0007)
   )
   
   (Partit 
 (Root (0.0 {parameter[44]}) (0.3 {parameter[45]}) (0.5 0.4) (0.8 0.16) (0.95 0.04) (1.0 0.00) (2.0 0.00))
 (Leaf (0.0 {[parameter[45] if parameter[45] + 0.28 < 1 else 1-0.28 for i in range(1)][0]}) (0.3 {[parameter[46] if parameter[46] + 0.42 < 1 else 1-0.42 for i in range(1)][0]}) (0.5 {[parameter[47] if parameter[47] + 0.47 < 1 else 1-0.47 for i in range(1)][0]}) (0.8 0.3) (0.95 0.25) (1.1 0.00) (2.0 0.00))
 (Stem (0.0 0.28) (0.3 0.42) (0.5 0.47) (0.8 0.25) (0.95 0.62) (1.0 0.54) (1.1 0.38) (1.5 0.00) (2.0 0.00))         
 (NNI_crit 1.40 []); RATJEN & KAGE: 1.38
 (NNI_inc 0.65 []); RATJEN & KAGE: 0.60
 (RSR (0.0 1.1) (0.05 1) (0.3 0.5) (1.0 0.5) (2.0 0.25))
   )
   
   
   (Prod (NCrop 0.4)
(E_Root     0.69)
(E_Leaf     {parameter[49]})
(E_Stem     0.66)
(E_SOrg     0.7)
(r_Root     0.015)
(r_Leaf     0.016)
(r_Stem     0.01)
(r_SOrg     {parameter[50]})
(ShldResC   0.4)
(ReMobilDS  1.1)
(ReMobilRt  0.03)
(Large_RtDR 0.05)
(ExfoliationFac 0.7)
(LfDR (0.0 0.00) (0.3 0.00) (0.5 0.015) (1.0 {parameter[51]}) (1.5 {parameter[52]}) (2.0 0.08))
(RtDR (0.0 0.00) (0.5 0.004) (1.0 0.02) (2.0 0.03))
)
   
   (CrpN 
   (PtRootCnc (0.2 0.025) (0.5 0.026) (0.8 0.018) (1.20 0.015) (2.00 0.014))
   (CrRootCnc (0.2 0.017) (0.5 0.019) (0.8 0.014) (1.20 0.012) (2.00 0.009))
   (NfRootCnc (0.2 0.009) (1.20 0.003) (2.00 0.003))    
   (PtLeafCnc (0.2 0.051) (0.5 0.048) (0.8 0.044) (1.00 0.038) (2.00 0.021))
   (CrLeafCnc (0.2 0.039) (0.5 0.036) (1.00 0.025) (2.00 0.014))
   (NfLeafCnc (0.2 0.01) (0.5 0.01) (1.00 0.006) (2.00 0.005))
   (PtStemCnc (0.2 0.038) (0.5 0.035) (0.8 0.02) (1.00 0.015) (2.00 0.011))
   (CrStemCnc (0.2 0.023) (0.5 0.013) (0.8 0.008) (1.00 0.007) (2.00 0.006))
   (NfStemCnc (0.2 0.01) (1.00 0.002) (2.00 0.002))
   (PtSOrgCnc (1.00 0.017) (2.00 0.0175))
   (CrSOrgCnc (1.00 0.0125) (2.00 0.0125))
   (NfSOrgCnc (1.00 0.0005) (2.00 0.0005))
   (TLLeafEff (0 0.9) (2 0.9))
   (TLRootEff (0 0.1) (2 0.1))
   )
 
(Harvest (EconomicYield_W 0.8) (EconomicYield_N 0.94))
)'''
 
    ir_normal = irrigation_threshold_calculator(rootdepth = 130, nfk_drop = parameter[17])
    ir_wet = irrigation_threshold_calculator(rootdepth = 130, nfk_drop = parameter[35])
    ir_dry = irrigation_threshold_calculator(rootdepth = 130, nfk_drop = parameter[53])
    irrigation_normal = f'''(defaction irrigate_normal activity
 (wait (and (after_mm_dd 4 14)
(before_mm_dd 7 31)
;;nfk_threshold as {parameter[17]}%
(not (soil_water_content_above {ir_normal} [mm] (from 0 [cm]) (to -130 [cm]))))) 
 (irrigate_overhead 25 [mm/h])
 (wait_days 4))'''
    irrigation_wet = f'''(defaction irrigate_wet activity
 (wait (and (after_mm_dd 4 14)
(before_mm_dd 7 31)
;;nfk_threshold as {parameter[35]}%
(not (soil_water_content_above {ir_wet} [mm] (from 0 [cm]) (to -130 [cm]))))) 
 (irrigate_overhead 25 [mm/h])
 (wait_days 4))'''
    irrigation_dry = f'''(defaction irrigate_dry activity
   (wait (and (after_mm_dd 4 14)
  (before_mm_dd 7 31)
  ;;nfk_threshold as {parameter[53]}%
  (not (soil_water_content_above {ir_dry} [mm] (from 0 [cm]) (to -130 [cm]))))) 
   (irrigate_overhead 25 [mm/h])
   (wait_days 4))'''
   
    management = f'''\n(defaction WinterWheat_custom_2005_management activity
	;; Winter Wheat 2005 management.
	(wait (at 2005 10 10))
	(sow "WinterWheat_normal"(seed 187.5 [kg w.w./ha]))
	(wait (at 2006 4 6))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2006 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 64.0 [kg N/ha]) ))
	(wait (at 2006 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2006 8 3))
	(harvest "WinterWheat_normal"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2005_management_auto activity
                (while (WinterWheat_custom_2005_management)
                	(repeat irrigate_normal)))

(defaction WinterWheat_custom_2006_management activity
	;; Winter Wheat 2006 management.
	(wait (at 2006 10 10))
	(sow "WinterWheat_wet"(seed 187.5 [kg w.w./ha]))
	(wait (at 2007 4 3))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2007 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 64.0 [kg N/ha]) ))
	(wait (at 2007 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2007 7 18))
	(harvest "WinterWheat_wet"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2006_management_auto activity
                (while (WinterWheat_custom_2006_management)
                	(repeat irrigate_wet)))

(defaction WinterWheat_custom_2007_management activity
	;; Winter Wheat 2007 management.
	(wait (at 2007 10 6))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2008 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2008 4 16))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2008 5 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2008 7 26))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2007_management_auto activity
                (while (WinterWheat_custom_2007_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2008_management activity
	;; Winter Wheat 2008 management.
	(wait (at 2008 10 20))
	(sow "WinterWheat_normal"(seed 187.5 [kg w.w./ha]))
	(wait (at 2009 3 6))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2009 4 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2009 5 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2009 8 4))
	(harvest "WinterWheat_normal"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2008_management_auto activity
                (while (WinterWheat_custom_2008_management)
                	(repeat irrigate_normal)))

(defaction WinterWheat_custom_2009_management activity
	;; Winter Wheat 2009 management.
	(wait (at 2009 10 15))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2010 3 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2010 4 29))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2010 6 1))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2010 8 21))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2009_management_auto activity
                (while (WinterWheat_custom_2009_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2010_management activity
	;; Winter Wheat 2010 management.
	(wait (at 2010 10 26))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2011 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2011 4 11))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2011 5 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2011 8 17))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2010_management_auto activity
                (while (WinterWheat_custom_2010_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2012_management activity
	;; Winter Wheat 2012 management.
	(wait (at 2012 10 18))
	(sow "WinterWheat_wet"(seed 187.5 [kg w.w./ha]))
	(wait (at 2013 3 5))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2013 4 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2013 5 24))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2013 6 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2013 8 3))
	(harvest "WinterWheat_wet"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2012_management_auto activity
                (while (WinterWheat_custom_2012_management)
                	(repeat irrigate_wet)))

(defaction WinterWheat_custom_2013_management activity
	;; Winter Wheat 2013 management.
	(wait (at 2013 10 2))
	(sow "WinterWheat_normal"(seed 187.5 [kg w.w./ha]))
	(wait (at 2014 2 25))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2014 3 31))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2014 4 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2014 5 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2014 8 5))
	(harvest "WinterWheat_normal"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2013_management_auto activity
                (while (WinterWheat_custom_2013_management)
                	(repeat irrigate_normal)))

(defaction WinterWheat_custom_2014_management activity
	;; Winter Wheat 2014 management.
	(wait (at 2014 10 7))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2015 3 9))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2015 4 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2015 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2015 5 25))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2015 6 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2015 8 7))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2014_management_auto activity
                (while (WinterWheat_custom_2014_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2015_management activity
	;; Winter Wheat 2015 management.
	(wait (at 2015 10 21))
	(sow "WinterWheat_wet"(seed 187.5 [kg w.w./ha]))
	(wait (at 2016 3 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2016 4 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2016 5 11))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2016 6 1))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2016 8 8))
	(harvest "WinterWheat_wet"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2015_management_auto activity
                (while (WinterWheat_custom_2015_management)
                	(repeat irrigate_wet)))

(defaction WinterWheat_custom_2016_management activity
	;; Winter Wheat 2016 management.
	(wait (at 2016 10 27))
	(sow "WinterWheat_wet"(seed 187.5 [kg w.w./ha]))
	(wait (at 2017 3 9))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2017 4 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2017 5 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2017 6 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2017 8 8))
	(harvest "WinterWheat_wet"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2016_management_auto activity
                (while (WinterWheat_custom_2016_management)
                	(repeat irrigate_wet)))

(defaction WinterWheat_custom_2017_management activity
	;; Winter Wheat 2017 management.
	(wait (at 2017 10 19))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2018 3 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2018 4 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 72.0 [kg N/ha]) ))
	(wait (at 2018 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2018 6 5))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2018 7 23))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2017_management_auto activity
                (while (WinterWheat_custom_2017_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2018_management activity
	;; Winter Wheat 2018 management.
	(wait (at 2018 10 19))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2019 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2019 4 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2019 5 9))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2019 5 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2019 7 23))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2018_management_auto activity
                (while (WinterWheat_custom_2018_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2019_management activity
	;; Winter Wheat 2019 management.
	(wait (at 2019 10 15))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2020 3 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2020 4 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2020 5 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2020 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2020 7 24))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2019_management_auto activity
                (while (WinterWheat_custom_2019_management)
                	(repeat irrigate_dry)))

(defaction WinterWheat_custom_2020_management activity
	;; Winter Wheat 2020 management.
	(wait (at 2020 10 15))
	(sow "WinterWheat_wet"(seed 187.5 [kg w.w./ha]))
	(wait (at 2021 3 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2021 4 13))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 65.0 [kg N/ha]) ))
	(wait (at 2021 5 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 87.0 [kg N/ha]) ))
	(wait (at 2021 7 30))
	(harvest "WinterWheat_wet"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2020_management_auto activity
                (while (WinterWheat_custom_2020_management)
                	(repeat irrigate_wet)))

(defaction WinterWheat_custom_2021_management activity
	;; Winter Wheat 2021 management.
	(wait (at 2021 10 28))
	(sow "WinterWheat_dry"(seed 187.5 [kg w.w./ha]))
	(wait (at 2022 3 14))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2022 4 13))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2022 5 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 62.0 [kg N/ha]) ))
	(wait (at 2022 7 26))
	(harvest "WinterWheat_dry"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2021_management_auto activity
                (while (WinterWheat_custom_2021_management)
                	(repeat irrigate_dry)))

;; Combine it.
        (defprogram Hamerstorf Daisy
            ;;Weather data
            (weather default "Hamerstorf.dwf")
            ;;Field to use
            (column Hamerstorf)'''
    management_years_calibration = f'''
;Simulation start and stop dates
    (time 2013 9 1)
    (stop 2022 7 31)
;Management start
        (manager activity
	;WinterWheat_custom_2005_management_auto
	;WinterWheat_custom_2006_management_auto
	;WinterWheat_custom_2007_management_auto
	;WinterWheat_custom_2008_management_auto
	;WinterWheat_custom_2009_management_auto
	;WinterWheat_custom_2010_management_auto
	;WinterWheat_custom_2012_management_auto
	WinterWheat_custom_2013_management_auto
	WinterWheat_custom_2014_management_auto
	WinterWheat_custom_2015_management_auto
	WinterWheat_custom_2016_management_auto
	WinterWheat_custom_2017_management_auto
	WinterWheat_custom_2018_management_auto
	WinterWheat_custom_2019_management_auto
	WinterWheat_custom_2020_management_auto
	WinterWheat_custom_2021_management_auto)
	 
;;Create these log files 
    	(output harvest
    		("Field water" (when monthly) (where "field_water.dlf"))
    		("Field water" (when monthly) (where "field_water.dlf"))
    		;("Field water" (when daily) (where "field_waterdaily.dlf"))
    		("Soil water" (when daily)
    			(where "soil_water30cm.dlf")
    			(from 0 [cm]) (to -30 [cm]))
    		;;("Soil water" (when daily)
    			;;(where "soil_water60cm.dlf")
    			;;(from -30 [cm]) (to -60 [cm]))
    		("Crop Production")))
;;Use it
    (run Hamerstorf)
     '''
    management_years_validation = f'''
;Simulation start and stop dates
    (time 2001 1 1)
    (stop 2022 7 31)
;Management start
        (manager activity
	;WinterWheat_custom_2005_management_auto
	;WinterWheat_custom_2006_management_auto
	WinterWheat_custom_2007_management_auto
	WinterWheat_custom_2008_management_auto
	WinterWheat_custom_2009_management_auto
	WinterWheat_custom_2010_management_auto
	WinterWheat_custom_2012_management_auto
	WinterWheat_custom_2013_management_auto
	WinterWheat_custom_2014_management_auto
	WinterWheat_custom_2015_management_auto
	WinterWheat_custom_2016_management_auto
	WinterWheat_custom_2017_management_auto
	WinterWheat_custom_2018_management_auto
	WinterWheat_custom_2019_management_auto
	WinterWheat_custom_2020_management_auto
	WinterWheat_custom_2021_management_auto)
	 
;;Create these log files 
    	(output harvest
    		("Field water" (when monthly) (where "field_water.dlf"))
    		("Field water" (when monthly) (where "field_water.dlf"))
    		;("Field water" (when daily) (where "field_waterdaily.dlf"))
    		("Soil water" (when daily)
    			(where "soil_water30cm.dlf")
    			(from 0 [cm]) (to -30 [cm]))
    		;;("Soil water" (when daily)
    			;;(where "soil_water60cm.dlf")
    			;;(from -30 [cm]) (to -60 [cm]))
    		("Crop Production")))
;;Use it
    (run Hamerstorf)
     '''
     
    
    if calibration:
        all_stmts = [initialize, soil_horizons, soil_column, crop_normal, crop_wet, crop_dry, 
                     irrigation_normal, irrigation_wet,irrigation_dry, management, management_years_calibration]
    else:
        all_stmts = [initialize, soil_horizons, soil_column, crop_normal, crop_wet, crop_dry, 
                     irrigation_normal, irrigation_wet,irrigation_dry, management, management_years_validation]
    run_txt = os.path.join(dir_, 'wwheat_nsgaII.dai')
    with open(run_txt, 'w') as exec_file:
        for text in all_stmts:
            exec_file.writelines(text)
            exec_file.writelines('\n')
            
def daisy_hame(working_dir,         #main directory
              exec_count,           #will go to directory where daisy will execute
              parameter_names,      #names of parameters
              parameter_block,      #2d parameter array over which daisy will loop through
              savedir_multiplier,   #keeps track of where to save each daisy dlf after execution
              calibration):         #true for calibration, false for validation
    
    #sends daisy to a dir where it will execute
    
    e_dir = os.path.join(working_dir, 'exec_dir',str(exec_count))
    s_dir_first_path = exec_count * savedir_multiplier
    count = 0
    for parameter_vector_count in range(parameter_block.shape[0]):
        parameter_vector = parameter_block[parameter_vector_count,:]
        print(parameter_vector)
        #print(parameter_vector, parameter_vector.shape)
        s_dir_path = s_dir_first_path + count
        count += 1
        s_dir = os.path.join(working_dir, 'save_dir', str(s_dir_path))
        #save the parameter vector used
        df = pd.DataFrame(parameter_vector.reshape(1,parameter_block.shape[1]), columns = parameter_names)
        parameter_txt_save_path = os.path.join(e_dir, 'parameters.txt')
        df.to_csv(parameter_txt_save_path, sep = ';', header = True, index = False)
        
        #create management file
        write_management_nsga_II(e_dir, parameter_vector, calibration)
        
        #executes daisy
        os.chdir(e_dir)
        try:
            print('Daisy Running')
            subprocess.run('daisy.exe wwheat_nsgaII.dai',
                shell=True, check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
        except:
            print('Daisy execution failed.')
        
        #copy simulation outputs and parameters used to savefile
        sim_outputs = ['wwheat_nsgaII.dai', 'parameters.txt','daisy.log',        #debugging if params wrong
                       'harvest.dlf', 'field_water.dlf', 
                       'soil_water30cm.dlf', 'soil_water60cm.dlf']
        for output in sim_outputs:
            source_path = os.path.join(e_dir, output)
            destination_path = os.path.join(s_dir, output)
            try:
                shutil.copy(source_path, destination_path)
            except:
                print(f'Parameterization problem')
                
def read_sm(sm_path, soil_depth,i, dates):
    
    df = pd.read_csv(sm_path, skiprows=20, header=[0, 1], sep='\t')
    df.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df.columns.values]
    df = df[df['mday'] != 0]
    df['Date'] = pd.to_datetime(df.rename(columns={'mday': 'day'})[['year', 'month', 'day']])
    df = df[df['Date'].isin(dates)]
    # print(df)
    df = df.loc[:, ['Date', 'Matrix water mm']]
    df[f'sim_sm_perc_{i}'] = df['Matrix water mm'] / soil_depth
    # df = df.loc[(df['Date'] > pd.to_datetime(start_date) 
    #              & (df['Date'] < pd.to_datetime(end_date))]
    
    # df = df.loc[:,['Date', f'sim_sm_perc_{i}']]
    sm = df[f'sim_sm_perc_{i}'].to_numpy()
    return sm

def read_harvest(harvest_dlf, years):
    df_harvest = pd.read_csv(harvest_dlf, skiprows = 9, header = [0,1], sep = '\t')
    df_harvest.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_harvest.columns.values]
    df_harvest = df_harvest[df_harvest['year'].isin(years)]
    harvest = df_harvest['sorg_DM Mg DM/ha'].to_numpy()
    return harvest


def read_irrigation(irrigation_dlf, years):
    df_irrigate = pd.read_csv(irrigation_dlf, skiprows = 26, header = [0,1], sep = '\t')
    df_irrigate.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_irrigate.columns.values]
    df_irrigate = df_irrigate[['year', 'month' , 'mday','Irrigation mm']]
    #df_irrigate = df_irrigate.iloc[5:].reset_index(drop= True)
    df_irrigate['crop_year'] = np.where(df_irrigate['month'] < 10, df_irrigate['year'],df_irrigate['year']+1)
    df_irrigate = df_irrigate[df_irrigate['crop_year'].isin(years)]
    df_irrigate = df_irrigate.groupby('crop_year')['Irrigation mm'].sum(min_count=1).reset_index()
    irrigation = df_irrigate['Irrigation mm'].to_numpy()
    return irrigation


    

def daisy_read_n_eval(main_dir, yld_path, sm_path, years, NSGA = True):
    #read observed values to calculate rmse
    #hardcoded index for calibration period
    
    
    #read observed data
    df_obsyld = pd.read_excel(yld_path, header = [0], sheet_name='ww_optimal')
    df_obsyld = df_obsyld [df_obsyld ['years'].isin(years)]
    obsyld = df_obsyld.iloc[:,1].to_numpy()                                  
    obsirr = df_obsyld.iloc[:,2].to_numpy()                                     
    
    df_obssm = pd.read_excel(sm_path, header = [0])
    df_obs30 = df_obssm.loc[:,['Date', 'obs moist_30 cm', 'year']]
    df_obs30 = df_obssm[df_obs30['year'].isin(years)]
    
    obssm = df_obs30['obs moist_30 cm'].to_numpy()
    sm_dates = df_obs30['Date'].to_numpy()
    
    #create list to store data
    
    yield_outputs = []
    irrigation_outputs = []
    sm_outputs = []
    parameters = []

    
    # print(sm_dates)
    # print(sm_dates.type())
    
    #read individual simulation
    for i in range(0,len(os.listdir(os.path.join(main_dir, 'save_dir')))):
        
        #creates file path of output to read
        parameter_txt = os.path.join(main_dir, 'save_dir', str(i),'parameters.txt')
        harvest_dlf = os.path.join(main_dir, 'save_dir', str(i),'harvest.dlf')
        irrigation_dlf = os.path.join(main_dir, 'save_dir', str(i),'field_water.dlf')
        sm30_dlf = os.path.join(main_dir, 'save_dir', str(i),'soil_water30cm.dlf')
        
        #parameter
        df_parameter = pd.read_csv(parameter_txt, header = [0], sep = ';')
        if i == 0:
            parameter_names = df_parameter.columns.to_list()
            print(parameter_names)
        parameters.append(df_parameter.to_numpy()[0,:])
        
        #harvest
        yield_this_model_run = read_harvest(harvest_dlf, years)
        yield_outputs.append(yield_this_model_run)
        
        #irrigate
        irrigation_this_model_run = read_irrigation(irrigation_dlf, years)
        irrigation_outputs.append(irrigation_this_model_run)
        
        #soil_moisture
        sm_30_this_model_run = read_sm(sm30_dlf, 300,i, sm_dates)
        sm_outputs.append(sm_30_this_model_run)
        #     #print(df_irrigate)
    
    
    #compile all simulations
    #parameter_array = np.array(parameters)
    yield_array = np.array(yield_outputs)
    irrigation_array = np.array(irrigation_outputs)
    sm_array = np.array(sm_outputs)
    print(yield_array.shape)
    
    #calculate_rmse here                                                        #define axis properly, #return in format suitable for pymoo
    rmse_yld = np.sqrt(np.mean((obsyld-yield_array)**2, axis = 1))
    rmse_irr = np.sqrt(np.mean((obsirr-irrigation_array)**2, axis = 1))
    rmse_sm =  np.sqrt(np.mean((obssm -sm_array)**2, axis = 1))
    
    # irrigation_iwue = np.where(irrigation_array == 0, 1, irrigation_array )
    # iwue = np.mean(yield_array/irrigation_iwue, axis = 1)
    F = np.column_stack((rmse_yld, rmse_irr, rmse_sm))
    
    if NSGA:
        # return F,G
        return F
    else:
        return yield_array, irrigation_array, sm_array

def daisy_kame(main_dir, param_names, param_blocks, yld_path, sm_path, years, NSGA, calibration):  
#this functions parameter blocks to function where a for loop               #
        #is run for all the parameters in the block. 
    #the parameters used are saved, and outputs of the daisy run are pasted
        #to save directory
    process_count = len(param_blocks)
    #print('Number of parameter blocks sent: ', process_count)
    save_dir_multiplier = param_blocks[0].shape[0]
    print('Saving files separated by an AP, a= 0 d = ', save_dir_multiplier)
    exec_dir = os.path.join(main_dir, 'exec_dir')                               #exec_dir = execute daisy here
    if not os.path.exists(exec_dir):
        os.makedirs(exec_dir)
    for i in range(process_count):
        sub_exec_dir = os.path.join(exec_dir, str(i))
        if not os.path.exists(sub_exec_dir):
            os.makedirs(sub_exec_dir)
    
    save_dir =  os.path.join(main_dir, 'save_dir')                              #create save dir here
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
        os.makedirs(save_dir)
    elif not os.path.exists(save_dir):
        os.makedirs(save_dir)
    total_save_dirs = sum([block.shape[0] for block in param_blocks])
    for j in range(total_save_dirs):
        sub_save_dir = os.path.join(save_dir, str(j))
        if not os.path.exists(sub_save_dir):
            os.makedirs(sub_save_dir)
        
    
    
    mc_kwrds = [(main_dir, k, param_names, param_blocks[k],save_dir_multiplier, calibration) 
                for k in range(process_count)]
    #print(mc_kwrds[0])
    with Pool(processes = process_count) as pool:
        try:
            results = pool.starmap(daisy_hame, mc_kwrds)
            for i, output in enumerate(results):
                    print(f"[{i}] Finished successfully:\n{output}", flush=True)
        except Exception as e:
            print(f"Error in multiprocessing: {e}", flush=True)                 #dasiy will complete execution and saving until this point
    F = daisy_read_n_eval(main_dir, yld_path, sm_path, years, NSGA)
    return F                          #main_dir= change with crop, param_blocks = changeable
    #this functions parameter blocks to function where a for loop               #
        #is run for all the parameters in the block. 
    #the parameters used are saved, and outputs of the daisy run are pasted
        #to save directory
    process_count = len(param_blocks)
    #print('Number of parameter blocks sent: ', process_count)
    save_dir_multiplier = param_blocks[0].shape[0]
    print('Saving files separated by an AP, a= 0 d = ', save_dir_multiplier)
    exec_dir = os.path.join(main_dir, 'exec_dir')                               #exec_dir = execute daisy here
    if not os.path.exists(exec_dir):
        os.makedirs(exec_dir)
    for i in range(process_count):
        sub_exec_dir = os.path.join(exec_dir, str(i))
        if not os.path.exists(sub_exec_dir):
            os.makedirs(sub_exec_dir)
    
    save_dir =  os.path.join(main_dir, 'save_dir')                              #create save dir here
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
        os.makedirs(save_dir)
    elif not os.path.exists(save_dir):
        os.makedirs(save_dir)
    total_save_dirs = sum([block.shape[0] for block in param_blocks])
    for j in range(total_save_dirs):
        sub_save_dir = os.path.join(save_dir, str(j))
        if not os.path.exists(sub_save_dir):
            os.makedirs(sub_save_dir)
        
    
    
    mc_kwrds = [(main_dir, k, param_names, param_blocks[k],save_dir_multiplier) 
                for k in range(process_count)]
    #print(mc_kwrds[0])
    with Pool(processes = process_count) as pool:
        try:
            results = pool.starmap(daisy_hame, mc_kwrds)
            for i, output in enumerate(results):
                    print(f"[{i}] Finished successfully:\n{output}", flush=True)
        except Exception as e:
            print(f"Error in multiprocessing: {e}", flush=True)                 #dasiy will complete execution and saving until this point
    F = daisy_read_n_eval(main_dir)
    return F


# if __name__ == "__main__":
    
#     #execute_pareto
#     nsgaII_iteration_save_path =r'C:\Users\Adhikari\Desktop\wwheat_stepwise\2_a_wwheat_nsga\pc\0'
    
#     pareto_parameters = pd.read_excel(os.path.join(nsgaII_iteration_save_path, '0_pareto_front.xlsx'))
#     pareto_parameters = pareto_parameters.iloc[:,:-3] #beware of +1 indexing error here  
#     p_names = list(pareto_parameters.columns)
#     print(p_names)
#     parameter_blocks = param_blocks(pareto_parameters.to_numpy())
    
#     #execute pareto front
#     # daisy_kame(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\2_a_wwheat_nsga\pc\0',
#     #             p_names, parameter_blocks)
    
#     yield_par, irr_par, sm_par = daisy_read_n_eval(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\2_a_wwheat_nsga\pc\0', NSGA = False)
#     years = np.array([y for y in range(2008,2020) if y!= 2012])
#     df_obsyld = pd.read_excel(r'C:\Users\Adhikari\Desktop\Thesis\Thesis Shared Cloud\all_crop_field_results.xlsx', header = [0])
#     df_obsyld = df_obsyld [df_obsyld ['years'].isin(years)]
#     obsyld = df_obsyld.iloc[:,1].to_numpy() 
#     obsirr = df_obsyld.iloc[:,2].to_numpy()
    
#     for pareto_ind in range(yield_par.shape[0]):
#         fig, ax = plt.subplots(figsize=(12, 8))
#         ax_irr = ax.twinx()
#         yield_ind = yield_par[pareto_ind,:]
#         irr_ind = irr_par[pareto_ind,:]
#         yield_min = np.max(yield_par, axis = 0)
#         scatter1 = ax.plot(years,yield_ind, label='Simulated Yield', marker='o', linestyle='-', zorder=3)
#         scatter4 = ax.plot(years, yield_ind/irr_ind *100/2,label='Simulated IWUE', marker='x', linestyle='-.', zorder=3)
#         scatter5 = ax.plot(years, obsyld/obsirr *100/2,label='Observed IWUE', marker='x', linestyle='-.', zorder=3)
#         bar1 = ax_irr.bar(years-0.125, irr_ind/10, width=0.25, alpha=0.5, label='Simulated Irrigation', zorder=2, color = 'orange')
#         scatter2 = ax.plot(years,obsyld, label='Observed Yield', marker='o', linestyle='-',zorder=3)
#         bar2 = ax_irr.bar(years+0.125, obsirr/10, width=0.25, alpha=0.5, label='Observed Irrigation', zorder=2, color = 'blue')
#         ax.set_xticks(years)
#         ax.set_ylabel('Yield [ton/ha]')
#         ax.set_yticks(np.arange(0,21,1))
#         ax_irr.set_ylabel('Irrigation [cm]')
#         ax_irr.set_yticks(np.arange(0,29,2.5))
#         plt.title(f'Pareto Parameter Set {pareto_ind}')
#         yld_rmse1 = np.sqrt(np.mean(yield_ind-obsyld)**2)
#         irr_rmse1 = np.sqrt(np.mean(irr_ind-obsirr)**2)
#         txt = f'yield rmse = {yld_rmse1}\nirrigation rmse = {irr_rmse1}'
#         #ax.text((2013, 10), txt, zorder=3)
#         lines1, labels1 = ax.get_legend_handles_labels()
#         lines2, labels2 = ax_irr.get_legend_handles_labels()
#         ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
#         print(yield_ind/irr_ind *100)
#         file = f'Pareto Parameter Set {pareto_ind}.png'
#         file_path = os.path.join(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\2_a_wwheat_nsga\pc\0', file)
#         plt.savefig(file_path)
#         plt.show()
        
    