import math as math
import numpy as np
import globals
#from globals import coef_from_txt
#from globals import coef_readed

def force_coef_body(mach,alfa,beta):

############################
    if 1:
# Costello-Sahu AIAA 2007-6582
#Spark data
        #Cx0 = 0.221
        #Cx2 = 5.0
        #Cna = 5.83
        #Cypa = 0
 #Estimation CFD/RBD
        Cx0 = 0.2387
        Cx2 = 5.933
        Cna = 5.644
        Cypa = 0
##########################
#Stahl-Costello-Sahu AIAA 2009-5715
#Synthetic
    if 0:
        Cx0 = 0.23
        Cx2 = 5
        Cna = 5.9
        Cypa = 0
    #Spark data
    #Cx0 = 0.221
    #Cx2 = 5.0
    #Cna = 5.83
    if globals.Forces_coef_from_txt:
        if not globals.Forces_coef_readed:
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_baranwonski.txt', delimiter=',', skiprows=3)
            data_raw = np.loadtxt('Resu_ref/body_coef_txt/bc_egipcio.txt', delimiter=',', skiprows=3)
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_NWU_104pg.txt', delimiter=',', skiprows=3)
            data = data_raw[0:-1]
            globals.data = data
            globals.data_raw = data_raw
            globals.Forces_coef_readed = True

        M = globals.data[:,0]
            
        # drag
        Cx0_exp = globals.data[:,1]
        Cx0 = np.interp(mach, M, Cx0_exp)

        Cx2_exp = globals.data[:,2]
        Cx2 = np.interp(mach, M, Cx2_exp)

        alfa_total2 = ((math.sin(beta))**2 + (math.cos(beta))**2*(math.sin(alfa))**2)
        Cx = Cx0 + Cx2*alfa_total2

        # normal
        Cna_exp = globals.data[:,3]
        Cna = np.interp(mach, M, Cna_exp)

        # magnus
        Cypa_exp = globals.data[:,4]
        Cypa = np.interp(mach, M, Cypa_exp)

    return [Cx, Cx2, Cna, Cypa]
