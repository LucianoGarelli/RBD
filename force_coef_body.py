import math as math
import numpy as np
import globals
#from globals import coef_from_txt
#from globals import coef_readed

def force_coef_body(mach,alfa,beta):

    #Spark data
    #Cx0 = 0.221
    #Cx2 = 5.0
    #Cna = 5.83
    if globals.Forces_coef_from_txt:
        if not globals.Forces_coef_readed:
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_baranwonski.txt', delimiter=',', skiprows=3)
            data = np.loadtxt('Resu_ref/body_coef_txt/bc_egipcio.txt', delimiter=',', skiprows=3)
            #data = np.loadtxt('Resu_ref/body_coef_txt/bc_NWU_104pg.txt', delimiter=',', skiprows=3)

            globals.data = data

            # Toma el primer valor del vector
            M_0 = data[0,0]
            Cx0_0 = data[0,1]
            Cx2_0 = data[0,2]
            Cna_0 = data[0,3]
            Cypa_0 = data[0,4]

            M =  M_0 
            Cx0 = Cx0_0 
            Cx2 = Cx2_0 
            Cna = Cna_0 
            Cypa = Cypa_0
            
            globals.Forces_coef_readed = True
        else:
            #
            M = globals.data[:,0]
            
            # drag
            Cx0_exp = globals.data[:,1]
            Cx0 = np.interp(mach, M, Cx0_exp)

            Cx2_exp = globals.data[:,2]
            Cx2 = np.interp(mach, M, Cx2_exp)

            #alfa_total2 = ((math.sin(beta))**2 + (math.cos(beta))**2*(math.sin(alfa))**2)
            #Cx = Cx0 + Cxd2*alfa_total2

            # normal
            Cna_exp = globals.data[:,3]
            Cna = np.interp(mach, M, Cna_exp)

            # magnus
            Cypa_exp = globals.data[:,4]
            Cypa = np.interp(mach, M, Cypa_exp)

        return [Cx0, Cx2, Cna, Cypa]
