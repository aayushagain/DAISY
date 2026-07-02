# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 16:53:39 2025

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


def param_blocks(param_values, mc_cores = None):
    if mc_cores == None:
        mc_cores = os.cpu_count()-2
    elif mc_cores > os.cpu_count():
        mc_cores = os.cpu_count()
    
    param_blocks = []
    mc_index = int(param_values.shape[0]/mc_cores)
    start = 0
    for i in range(mc_cores):
        if i < mc_cores-1:
            end = start + mc_index
            print(start,end)
            #param_small_index.append([start, end])
            param_blocks.append(param_values[start:end])
            start = end
        else:                                                                       #bypass if params cannot be divided into 16 equal parts
            end = param_values.shape[0]
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
        elif rootdepth > low:
            weight_list[i] = rootdepth - low
        else:
            weight_list[i] = 0
    weight = np.array(weight_list)
    # --- Threshold moisture calculation ---
    oall = pwp + nfk_drop/100 * (fc-pwp)
    wc_thres = (oall * weight).sum()*10
    return wc_thres

def write_management_morris_screening(dir_, parameter):
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
            (dry_bulk_density 1.64 [g/cm^3]))
    '''
    
    soil_column = f'''(defcolumn Hamerstorf default
          "Data collected by Hamerstorf."
          (Bioclimate default (pet PM))
          (Soil (horizons (-30 [cm] Ap) (-50 [cm] B) (-60 [cm] C) (-95 [cm] D) (-200 [cm] E))
                (MaxRootingDepth 100 [cm]))
          (Movement vertical(Geometry (zplus -2.5 -5 -10 -15 -20 -30 -40 -50 -60 -70 -80 -90  -100 -125 -150 -200 [cm])))
          (Groundwater deep))
    '''
    
    crop = f''';Creating crop parameters  
  (defcrop "WinterWheat_custom" "Winter Wheat"
    
   (Devel default 
     (EmrTSum {parameter[0]}) (DSRate1 {parameter[1]}) (DSRate2 {parameter[2]})
     (DS_Emr 0.01)
     (TempEff1 (-10.0  0.01) (0.00  0.01) (20.0  0.90) (25.0  1.00) (35.0  1.20))
     (PhotEff1 (10.0  0.29) (11.0  0.55) (12.0  0.75) (13.0  0.89) (14.0  1.00) (15.0  1.08) (16.0  1.14) (17.0  1.18) (24.0  1.18))
     (TempEff2 ( 0.0  0.00) (10.0  0.14) (15.0  0.66) (20.0  1.00) (25.0  1.23))
    )
    
   (Vernal default 
     (DSLim {parameter[3]}) (TaLim {parameter[4]}) (TaSum {parameter[5]})
    )
     
   (LeafPhot original 
     (Fm {parameter[6]}) 
     (Qeff   {parameter[7]})
     (TempEff (-20.0  0.00) (4.00  0.00) (10.0  1.00) (25.0  1.00) (35.0  0.01) (50.0  0.00))
    )
     
   (Canopy
     (SpLAI    {parameter[8]})
     (PARref   {parameter[9]})
     (PARext   {parameter[10]})
     (PARrel   {parameter[11]})
     (EPext    {parameter[12]})
     (LAIDist0 0.00 0.00 1.00)
     (LAIDist1 0.00 {parameter[13]} {parameter[14]})
     (HvsDS   (0.00    1)    (1.00  80)    (2.00 100))
    )
   
    (Root (DptEmr     {parameter[15]})
     (PenPar1    {parameter[16]})
     (PenPar2     {parameter[17]})
     (MaxPen      {parameter[18]})
     (Rad        {parameter[19]})
     (h_wp       -15000)
     (Rxylem     {parameter[20]})
     (MxNH4Up     {parameter[21]})
     (MxNO3Up     {parameter[22]})
    )
    
    (Partit 
      (Root (0.0 {parameter[23]}) (0.3 {parameter[24]}) (0.5 {parameter[25]}) (0.8 {parameter[26]}) (0.95 {parameter[27]}) (1.0 0.00) (2.0 0.00))
      (Leaf (0.0 {parameter[28]}) (0.3 {parameter[29]}) (0.5 {parameter[30]}) (0.8 {parameter[31]}) (0.95 {parameter[32]}) (1.1 0.00) (2.0 0.00))
      (Stem (0.0 {1-parameter[28]}) (0.3 {1-parameter[29]}) (0.5 {1-parameter[30]}) (0.8 {parameter[32]}) (0.95 {parameter[33]}) (1.0 {parameter[34]}) (1.1 {parameter[35]}) (1.5 0.00) (2.0 0.00))         
      (NNI_crit 1.40 []); RATJEN & KAGE: 1.38
      (NNI_inc 0.65 []); RATJEN & KAGE: 0.60
      (RSR (0.0 1.1) (0.05 1) (0.3 0.5) (1.0 {parameter[36]}) (2.0 {parameter[37]}))
    )
    
    
    (Prod (NCrop {parameter[38]})
         (E_Root     {parameter[39]})
         (E_Leaf     {parameter[40]})
         (E_Stem     {parameter[41]})
         (E_SOrg     {parameter[42]})
         (r_Root     {parameter[43]})
         (r_Leaf     {parameter[44]})
         (r_Stem     {parameter[45]})
         (r_SOrg     {parameter[46]})
         (ShldResC   {parameter[47]})
         (ReMobilDS  {parameter[48]})
         (ReMobilRt  {parameter[49]})
         (Large_RtDR {parameter[50]})
         (ExfoliationFac {parameter[51]})
         (LfDR (0.0 0.00) (0.3 0.00) (0.5 {parameter[52]}) (1.0 {parameter[53]}) (1.5 {parameter[54]}) (2.0 {parameter[55]}))
         (RtDR (0.0 0.00) (0.5 {parameter[56]}) (1.0 {parameter[57]}) (2.0 {parameter[58]}))
         )
    
    (CrpN 
        (PtRootCnc (0.2 {parameter[59]}) (0.5 {parameter[60]}) (0.8 {parameter[61]}) (1.20 {parameter[62]}) (2.00 {parameter[63]}))
        (CrRootCnc (0.2 {parameter[64]}) (0.5 {parameter[65]}) (0.8 {parameter[66]}) (1.20 {parameter[67]}) (2.00 {parameter[68]}))
        (NfRootCnc (0.2 {parameter[69]}) (1.20 {parameter[70]}) (2.00 {parameter[71]}))    
        (PtLeafCnc (0.2 {parameter[72]}) (0.5 {parameter[73]}) (0.8 {parameter[74]}) (1.00 {parameter[75]}) (2.00 {parameter[76]}))
        (CrLeafCnc (0.2 {parameter[77]}) (0.5 {parameter[78]}) (1.00 {parameter[79]}) (2.00 {parameter[80]}))
        (NfLeafCnc (0.2 {parameter[81]}) (0.5 {parameter[82]}) (1.00 {parameter[83]}) (2.00 {parameter[84]}))
        (PtStemCnc (0.2 {parameter[85]}) (0.5 {parameter[86]}) (0.8 {parameter[87]}) (1.00 {parameter[88]}) (2.00 {parameter[89]}))
        (CrStemCnc (0.2 {parameter[90]}) (0.5 {parameter[91]}) (0.8 {parameter[92]}) (1.00 {parameter[93]}) (2.00 {parameter[94]}))
        (NfStemCnc (0.2 {parameter[95]}) (1.00 {parameter[96]}) (2.00 {parameter[97]}))
        (PtSOrgCnc (1.00 {parameter[98]}) (2.00 {parameter[99]}))
        (CrSOrgCnc (1.00 {parameter[100]}) (2.00 {parameter[101]}))
        (NfSOrgCnc (1.00 {parameter[102]}) (2.00 {parameter[103]}))
        (TLLeafEff (0 {parameter[104]}) (2 {parameter[104]}))
        (TLRootEff (0 {parameter[105]}) (2 {parameter[105]}))
        )
  
  (Harvest (EconomicYield_W {parameter[106]}) (EconomicYield_N {parameter[107]}))
  )'''
    irrigation_threshold = irrigation_threshold_calculator(rootdepth = 80, nfk_drop = parameter[108]) 
   
    irrigation = f'''(defaction irrigate_auto activity
    (wait (and (after_mm_dd 4 14)
           (before_mm_dd 7 31)
           (not (soil_water_content_above {irrigation_threshold} [mm] (from 0 [cm]) (to -{80} [cm]))))) 
    (irrigate_overhead 25 [mm/h])
    (wait_days 4))'''
    
    management = f'''\n(defaction WinterWheat_custom_2005_management activity
	;; Winter Wheat 2005 management.
	(wait (at 2005 10 10))
	(sow "WinterWheat_custom")
	(wait (at 2006 4 6))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2006 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 64.0 [kg N/ha]) ))
	(wait (at 2006 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2006 8 3))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2005_management_auto activity
                (while (WinterWheat_custom_2005_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2006_management activity
	;; Winter Wheat 2006 management.
	(wait (at 2006 10 10))
	(sow "WinterWheat_custom")
	(wait (at 2007 4 3))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2007 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 64.0 [kg N/ha]) ))
	(wait (at 2007 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2007 7 18))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2006_management_auto activity
                (while (WinterWheat_custom_2006_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2007_management activity
	;; Winter Wheat 2007 management.
	(wait (at 2007 10 6))
	(sow "WinterWheat_custom")
	(wait (at 2008 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2008 4 16))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2008 5 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2008 7 26))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2007_management_auto activity
                (while (WinterWheat_custom_2007_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2008_management activity
	;; Winter Wheat 2008 management.
	(wait (at 2008 10 20))
	(sow "WinterWheat_custom")
	(wait (at 2009 3 6))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2009 4 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2009 5 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2009 8 4))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2008_management_auto activity
                (while (WinterWheat_custom_2008_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2009_management activity
	;; Winter Wheat 2009 management.
	(wait (at 2009 10 15))
	(sow "WinterWheat_custom")
	(wait (at 2010 3 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2010 4 29))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2010 6 1))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2010 8 21))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2009_management_auto activity
                (while (WinterWheat_custom_2009_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2010_management activity
	;; Winter Wheat 2010 management.
	(wait (at 2010 10 26))
	(sow "WinterWheat_custom")
	(wait (at 2011 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2011 4 11))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2011 5 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2011 8 17))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2010_management_auto activity
                (while (WinterWheat_custom_2010_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2012_management activity
	;; Winter Wheat 2012 management.
	(wait (at 2012 10 18))
	(sow "WinterWheat_custom")
	(wait (at 2013 3 5))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2013 4 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2013 5 24))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2013 6 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2013 8 3))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2012_management_auto activity
                (while (WinterWheat_custom_2012_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2013_management activity
	;; Winter Wheat 2013 management.
	(wait (at 2013 10 2))
	(sow "WinterWheat_custom")
	(wait (at 2014 2 25))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2014 3 31))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2014 4 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2014 5 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2014 8 5))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2013_management_auto activity
                (while (WinterWheat_custom_2013_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2014_management activity
	;; Winter Wheat 2014 management.
	(wait (at 2014 10 7))
	(sow "WinterWheat_custom")
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
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2014_management_auto activity
                (while (WinterWheat_custom_2014_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2015_management activity
	;; Winter Wheat 2015 management.
	(wait (at 2015 10 21))
	(sow "WinterWheat_custom")
	(wait (at 2016 3 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2016 4 20))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2016 5 11))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2016 6 1))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2016 8 8))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2015_management_auto activity
                (while (WinterWheat_custom_2015_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2016_management activity
	;; Winter Wheat 2016 management.
	(wait (at 2016 10 27))
	(sow "WinterWheat_custom")
	(wait (at 2017 3 9))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2017 4 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2017 5 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2017 6 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2017 8 8))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2016_management_auto activity
                (while (WinterWheat_custom_2016_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2017_management activity
	;; Winter Wheat 2017 management.
	(wait (at 2017 10 19))
	(sow "WinterWheat_custom")
	(wait (at 2018 3 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2018 4 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 72.0 [kg N/ha]) ))
	(wait (at 2018 5 8))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2018 6 5))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2018 7 23))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2017_management_auto activity
                (while (WinterWheat_custom_2017_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2018_management activity
	;; Winter Wheat 2018 management.
	(wait (at 2018 10 19))
	(sow "WinterWheat_custom")
	(wait (at 2019 3 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2019 4 4))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2019 5 9))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 30.0 [kg N/ha]) ))
	(wait (at 2019 5 22))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2019 7 23))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2018_management_auto activity
                (while (WinterWheat_custom_2018_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2019_management activity
	;; Winter Wheat 2019 management.
	(wait (at 2019 10 15))
	(sow "WinterWheat_custom")
	(wait (at 2020 3 18))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2020 4 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 60.0 [kg N/ha]) ))
	(wait (at 2020 5 7))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 50.0 [kg N/ha]) ))
	(wait (at 2020 6 2))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 40.0 [kg N/ha]) ))
	(wait (at 2020 7 24))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2019_management_auto activity
                (while (WinterWheat_custom_2019_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2020_management activity
	;; Winter Wheat 2020 management.
	(wait (at 2020 10 15))
	(sow "WinterWheat_custom")
	(wait (at 2021 3 10))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 70.0 [kg N/ha]) ))
	(wait (at 2021 4 13))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 65.0 [kg N/ha]) ))
	(wait (at 2021 5 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 87.0 [kg N/ha]) ))
	(wait (at 2021 7 30))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2020_management_auto activity
                (while (WinterWheat_custom_2020_management)
                	(repeat irrigate_auto)))

(defaction WinterWheat_custom_2021_management activity
	;; Winter Wheat 2021 management.
	(wait (at 2021 10 28))
	(sow "WinterWheat_custom")
	(wait (at 2022 3 14))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2022 4 13))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 80.0 [kg N/ha]) ))
	(wait (at 2022 5 21))
	(fertilize (AmmoniumNitrate (NH4_fraction 0.750 []) (weight 62.0 [kg N/ha]) ))
	(wait (at 2022 7 26))
	(harvest "WinterWheat_custom"))

;;Combining management and auto irrigation 
(defaction WinterWheat_custom_2021_management_auto activity
                (while (WinterWheat_custom_2021_management)
                	(repeat irrigate_auto)))

;; Combine it.
        (defprogram Hamerstorf Daisy
            ;;Weather data
            (weather default "Hamerstorf.dwf")
            ;;Field to use
            (column Hamerstorf)

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
    		("Soil water" (when daily)
    			(where "soil_water30cm.dlf")
    			(from 0 [cm]) (to -30 [cm]))
    		("Crop Production")))
;;Use it
    (run Hamerstorf)
     '''
    all_stmts = [initialize, soil_horizons, soil_column, crop, irrigation, management]
    run_txt = os.path.join(dir_, 'daisy_morris.dai')
    with open(run_txt, 'w') as exec_file:
        for text in all_stmts:
            exec_file.writelines(text)
            exec_file.writelines('\n')
            
def daisy_hame(working_dir,         #main directory
              exec_count,           #will go to directory where daisy will execute
              parameter_names,      #names of parameters
              parameter_block,      #2d parameter array over which daisy will loop through
              savedir_multiplier):  #keeps track of where to save each daisy dlf after execution
    
    #sends daisy to a dir where it will execute
    
    e_dir = os.path.join(working_dir, 'exec_dir',str(exec_count))
    s_dir_first_path = exec_count * savedir_multiplier
    count = 0
    for parameter_vector_count in range(parameter_block.shape[0]):
        parameter_vector = parameter_block[parameter_vector_count,:]
        #print(parameter_vector, parameter_vector.shape)
        s_dir_path = s_dir_first_path + count
        count += 1
        s_dir = os.path.join(working_dir, 'save_dir', str(s_dir_path))
        #save the parameter vector used
        df = pd.DataFrame(parameter_vector.reshape(1,parameter_block.shape[1]), columns = parameter_names)
        parameter_txt_save_path = os.path.join(e_dir, 'parameters.txt')
        df.to_csv(parameter_txt_save_path, sep = ';', header = True, index = False)
        
        #create management file
        write_management_morris_screening(e_dir, parameter_vector)
        
        #executes daisy
        os.chdir(e_dir)
        try:
            print('Daisy Running')
            subprocess.run('daisy.exe daisy_morris.dai',
                shell=True, check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)
        except:
            print('Daisy execution failed.')
        
        #copy simulation outputs and parameters used to savefile
        sim_outputs = ['daisy_morris.dai', 'parameters.txt','daisy.log',        #debugging if params wrong
                       'harvest.dlf', 'field_water.dlf', 
                       'soil_water30cm.dlf', 'soil_water60cm.dlf']
        for output in sim_outputs:
            source_path = os.path.join(e_dir, output)
            destination_path = os.path.join(s_dir, output)
            try:
                shutil.copy(source_path, destination_path)
            except:
                print(f'Parameterization problem')
                        



def daisy_kame(main_dir, param_names, param_blocks):                            #main_dir= change with crop, param_blocks = changeable
    #this functions parameter blocks to function where a for loop 
        #is run for all the parameters in the block. 
    #the parameters used are saved, and outputs of the daisy run are pasted
        #to save directory
    process_count = len(param_blocks)
    save_dir_multiplier = param_blocks[0].shape[0]
    
    exec_dir = os.path.join(main_dir, 'exec_dir')                               #exec_dir = execute daisy here
    if not os.path.exists(exec_dir):
        os.makedirs(exec_dir)
    for i in range(process_count):
        sub_exec_dir = os.path.join(exec_dir, str(i))
        if not os.path.exists(sub_exec_dir):
            os.makedirs(sub_exec_dir)
    
    save_dir =  os.path.join(main_dir, 'save_dir')                              #create save dir here
    if not os.path.exists(save_dir ):
        os.makedirs(save_dir )    
    total_save_dirs = sum([block.shape[0] for block in param_blocks])
    for j in range(total_save_dirs):
        sub_save_dir = os.path.join(save_dir, str(j))
        if not os.path.exists(sub_save_dir):
            os.makedirs(sub_save_dir)
        
    
    
    mc_kwrds = [(main_dir, k, param_names, param_blocks[k],save_dir_multiplier) 
                for k in range(process_count)]
    print(mc_kwrds)
    with Pool(processes = process_count) as pool:
        results = [pool.apply_async(daisy_hame, args=mc_kwrd) for mc_kwrd in mc_kwrds]
        
        for i,res in enumerate(results):
            print("In the results territory now",res, flush= True)
            try:
                output = res.get()
                print(f"[{i}] Finished successfully:\n{output}", flush=True)
            except:
                print(f"[{i}] Failed", flush=True)

def read_n_save_outputs(main_dir):
    for i in range(len(os.listdir(os.path.join(main_dir, 'save_dir')))):
        
        #reates file path of output to read
        parameter_txt = os.path.join(main_dir, 'save_dir', str(i),'parameters.txt')
        harvest_dlf = os.path.join(main_dir, 'save_dir', str(i),'harvest.dlf')
        irrigation_dlf = os.path.join(main_dir, 'save_dir', str(i),'field_water.dlf')
        #print(parameter_txt)
        
        #when reading the first output, create a df skeleton to which rest output will be concatenated to 
        if i == 0: 
            df_parameter = pd.read_csv(parameter_txt, header = [0], sep = ';')
            
            df_harvest = pd.read_csv(harvest_dlf, skiprows = 9, header = [0,1], sep = '\t')
            df_harvest.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_harvest.columns.values]
            df_yield = df_harvest[['year','sorg_DM Mg DM/ha']]
            df_yield = df_yield.rename(columns={'sorg_DM Mg DM/ha': i})
            
            df_irrigate = pd.read_csv(irrigation_dlf, skiprows = 26, header = [0,1], sep = '\t')
            df_irrigate.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_irrigate.columns.values]
            df_irrigate = df_irrigate[['year', 'month' , 'mday','Irrigation mm']]
            df_irrigate = df_irrigate.iloc[5:].reset_index(drop= True)
            df_irrigate['crop_year'] = np.where(df_irrigate['month'] <= 8, df_irrigate['year'],
                                            np.where(df_irrigate['month'] >= 10, df_irrigate['year']+1,df_irrigate['year']))
            df_irrigate = df_irrigate.groupby('crop_year')['Irrigation mm'].sum(min_count=1).reset_index()
            df_irrigate = df_irrigate.rename(columns={'Irrigation mm': i})
            print(df_irrigate)
        #reads the output and concatenates it to the first output read
        if i > 0:
            df = pd.read_csv(parameter_txt, header = [0], sep = ';')
            df_parameter = pd.concat([df_parameter, df])
            
            df_h = pd.read_csv(harvest_dlf, skiprows = 9, header = [0,1], sep = '\t')
            df_h.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_h.columns.values]
            # df_yield.loc[:,i] = df_h['sorg_DM Mg DM/ha'].values
            df_h = df_h.rename(columns={'sorg_DM Mg DM/ha': i})
            df_yield = pd.concat([df_yield, df_h[i]], axis = 1).reset_index(drop=True)
            # df_yield.merge(df_h[[str(i)]],left_index=True, right_index=True, how='left')
    
            df_i = pd.read_csv(irrigation_dlf, skiprows = 26,header = [0,1], sep = '\t')
            df_i.columns = [col[0] if 'Unnamed' in col[1] else ' '.join(col) for col in df_i.columns.values]
            df_i = df_i[['year', 'month' , 'mday','Irrigation mm']]
            df_i = df_i.iloc[5:].reset_index(drop= True)
            df_i['crop_year'] = np.where(df_i['month'] <= 8, df_i['year'],
                                            np.where(df_i['month'] >= 10, df_i['year']+1, df_i['year']))
            df_i = df_i.groupby('crop_year')['Irrigation mm'].sum(min_count=1).reset_index()
            df_i = df_i.rename(columns={'Irrigation mm': i})
            # print(df_i.head(10))
            df_irrigate = pd.concat([df_irrigate, df_i[i]], axis = 1).reset_index(drop=True)
            #print(df_irrigate)
        #remove the first 5 years of output if required
        

        #saves all read output as csv to a directory
        parameter_save_file_path = os.path.join(main_dir, 'parameters_combined.txt')
        df_parameter.to_csv(parameter_save_file_path, index = 0)

        yield_save_file_path = os.path.join(main_dir, 'yield.txt')
        df_yield_transposed = df_yield.T 
        df_yield_transposed.to_csv(yield_save_file_path, index = 0)
        
        irrigate_save_file_path = os.path.join(main_dir, 'irrigation.txt')
        df_irrigate_transposed = df_irrigate.T
        df_irrigate_transposed.to_csv(irrigate_save_file_path, index = 0)
    #return df_parameter, df_yield, df_yield_transposed, df_irrigate, df_irrigate_transposed         
    
if __name__ == "__main__":
    df = pd.read_excel(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\1_a_wwheat_morris\wwheat_morris.xlsx', 
                     header = [0])
    all_parameter = df.to_numpy()
    parameter_names = df.columns.tolist()
    parameter_blocks = param_blocks(all_parameter)
    daisy_kame(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\1_a_wwheat_morris', parameter_names , parameter_blocks)
    read_n_save_outputs(r'C:\Users\Adhikari\Desktop\wwheat_stepwise\1_a_wwheat_morris')
    
    
