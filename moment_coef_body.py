import math as math
import numpy as np
from scipy import interpolate
import globals
from procesar_magnus_moment_coef import procesar_magnus_moment_coef

def moment_coef_body(mach,alfa,beta):

############################
    if 0:
# Costello-Sahu AIAA 2007-6582
#Spark data
        #Cma = -12.6
        #Cmq = -196
        #Clp = -2.71
        #Cnpa = 0
#Estimation CFD/RBD
        Cma = -13.8278
        Cmq = -134.4
        Clp = -3.379
        Cnpa = 0.0
##########################
#Stahl-Costello-Sahu AIAA 2009-5715
#Synthetic
    if 0:
        Cma = -8.9
        Cmq = -198
        Clp = -2.6
        Cnpa = 0


    if globals.Moments_coef_from_txt:
        if not globals.Moments_coef_readed:
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_baranwonski.txt', delimiter=',', skiprows=3)
            data_raw = np.loadtxt('Resu_ref/body_coef_txt/bc_egipcio.txt', delimiter=',', skiprows=3)
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_NWU_104pg.txt', delimiter=',', skiprows=3)
            data = data_raw[0:-1]
            globals.data = data
            globals.data_raw = data_raw
            globals.Moments_coef_readed = True
            procesar_magnus_moment_coef()

        M = globals.data[:,0]
            
        # overturning also known as pitch yaw moment
        Cma_exp = globals.data[:,6]
        Cma = -1*np.interp(mach, M, Cma_exp)

        # rolling damping
        Clp_exp = globals.data[:,5]
        Clp = np.interp(mach, M, Clp_exp)
        # checkiar el 2

        # pitch yaw damping
        Cmq_exp = globals.data[:,7]
        Cmq = np.interp(mach, M, Cmq_exp)

        # magnus, tabla de doble entrada
        #
        # Como lo tratamos ?
        # incorporamos globals.Cnpa_proce en globals, hay que ver si lo dejamos o lo sacamos

        alfa_total2 = ((math.sin(beta))**2 + (math.cos(beta))**2*(math.sin(alfa))**2)
        alfa_total = np.sqrt(alfa_total2)
        Cnpa = interpolate.griddata((globals.Mpa,globals.ang),globals.Cnpa_proce,(mach,alfa_total),method='linear')  # va alfa_total ac'a ??
        #Cm_p_alfa = -1*interpolate.griddata((Mpa,alfa_mp),Cm_p_alfa_exp,(mach,alfa_total2),method='linear')



    return [Cma, Cmq, Clp, Cnpa]
