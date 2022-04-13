import math as math
import numpy as np
from scipy import interpolate
import globals
from procesar_magnus_moment_coef import procesar_magnus_moment_coef

def moment_coef(mach,alfa,beta):
    Clp=0
    Cm_alfa=0
    Cm_p_alfa=0
    Cm_q=0
    Cn_beta=0
    Cn_r=0
    ## Datos Ejemplo Carluccio pag. 244
    if 0 :
        ## Clp
        Clp = -0*0.01

        ## Cm_alfa
        Cm_alfa = 0*2.36

        ## Cm_p_alfa
        Cm_p_alfa = 0*0.02

        ## Cm_q + Cm_dot(alfa) -> Pitching damping moment debido a q + dot(alfa)
        Cm_q = -0*16.2



    #######################
    # Ejemplo .308" 168 grain Apendix A pag 217
    if 0:
        M = [0, 0.5, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]

        # Cl_p
        Clp_exp = [-0.015, -0.0125, -0.0108, -0.0107, -0.0105, -0.0103, -0.0100, -0.0099, -0.0098, -0.0095,
                   -0.0088, -0.0083, -0.008, -0.0075, -0.0073, -0.0068]

        Clp = 2*np.interp(mach, M, Clp_exp)

        ## Cm_alfa
        Cm_alfa0_exp = [3.05, 3.26, 3.38, 3.4, 3.43, 3.45, 3.24, 3.17, 3.15, 3.12, 3.06, 2.98, 2.88, 2.79, 2.69, 2.56]

        Cm_alfa0 = np.interp(mach, M, Cm_alfa0_exp)

        Md2 = [0, 0.95, 1.0, 1.05, 2.5]

        Cm_d2_exp = [-4.3, -4.3, -4.35, -4.4, -4.4]

        Cm_d2 = np.interp(mach, Md2, Cm_d2_exp)

        alfa_total2 = ((math.sin(beta)) ** 2 + (math.cos(beta)) ** 2 * (math.sin(alfa)) ** 2)
        Cm_alfa = Cm_alfa0 + Cm_d2 * alfa_total2

        # Cm_q + Cm_dot(alfa)

        Mq = [0, 1.05, 1.1, 1.2, 1.4, 1.6, 2.5]
        Cm_q_exp = [1.2, 1.2, 0.5, -3.6, -7.3, -8.2, -8.2]

        Cm_q = 2*np.interp(mach, Mq, Cm_q_exp)

        ## Cm_p_alfa
        Mpa = [0,0,0, 0.9,0.9,0.9, 1.1,1.1,1.1, 1.4,1.4,1.4, 1.7,1.7,1.7, 2.5,2.5,2.5]

        alfa_mp = [0,0.008895,0.12185, 0,0.008895,0.12185, 0,0.0056,0.12185, 0,0.003016,0.12185,
                   0,0.00171,0.12185, 0,0.00134,0.12185]

        Cm_p_alfa_exp = [-2.6, 0.06, 0.06, -2.6, 0.06, 0.06, -1.35, 0.05,0.05, -0.51,0.24,0.24,
                         -0.33,0.1,0.1, -0.33,0.1,0.1 ]

        #f = interpolate.interp2d(Mpa, alfa_mp, Cm_p_alfa_exp)
        Cm_p_alfa = -1*interpolate.griddata((Mpa,alfa_mp),Cm_p_alfa_exp,(mach,alfa_total2),method='linear')

    if globals.Moments_coef_from_txt:
        if not globals.Moments_coef_readed:
            # data = np.loadtxt('Resu_ref/body_coef_txt/bc_baranwonski.txt', delimiter=',', skiprows=3)
            #data_raw = np.loadtxt('Resu_ref/body_coef_txt/bc_egipcio.txt', delimiter=',', skiprows=3)
            data_raw = np.loadtxt('Resu_ref/body_coef_txt/Coef_estim_from_MHE_F01_unificated.txt', delimiter=',', skiprows=3)
            # data = np.loadtxt('Resu_ref/body_coef_txt/bc_NWU_104pg.txt', delimiter=',', skiprows=3)
            data = data_raw[0:-1]
            globals.data = data
            globals.data_raw = data_raw
            globals.Moments_coef_readed = True
            procesar_magnus_moment_coef()

        M = globals.data[:, 0]

        # overturning also known as pitch yaw moment
        Cma_exp = globals.data[:, 6]
        Cm_alfa = -1 * np.interp(mach, M, Cma_exp)

        # rolling damping
        Clp_exp = globals.data[:, 5]
        Clp = np.interp(mach, M, Clp_exp)/2
        # BRL = NACA/2

        # pitch yaw damping
        Cmq_exp = globals.data[:, 7]
        Cm_q = np.interp(mach, M, Cmq_exp)
        # BRL = NACA/2
        # magnus, tabla de doble entrada
        #
        # Como lo tratamos ?
        # incorporamos globals.Cnpa_proce en globals, hay que ver si lo dejamos o lo sacamos

        alfa_total2 = ((math.sin(beta)) ** 2 + (math.cos(beta)) ** 2 * (math.sin(alfa)) ** 2)
        alfa_total = np.sqrt(alfa_total2)
        #Cm_p_alfa = interpolate.griddata((globals.Mpa, globals.ang), globals.Cnpa_proce, (mach, alfa_total),
        #                            method='linear')
        Cm_p_alfa = globals.interp(mach, alfa_total)
        # BRL = NACA/2
        # Cm_p_alfa = -1*interpolate.griddata((Mpa,alfa_mp),Cm_p_alfa_exp,(mach,alfa_total2),method='linear')

        ## Debido a simetria de revolucion
        ## Cn_beta
        Cn_beta = -Cm_alfa

        ##Cn_r
        Cn_r = Cm_q

    return [Clp, Cm_alfa, Cm_p_alfa, Cm_q, Cn_beta, Cn_r]  # return[damping,]
