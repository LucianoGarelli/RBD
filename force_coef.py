import math as math
import numpy as np
import globals

def force_coef(mach, alfa, beta):

    Cd=0.0
    CL_alfa=0.0
    Cn_p_alfa=0
    Cn_q_alfa=0

    #######################
    # Datos Ejemplo Carlucci pag. 244
    # CD Coeficiente de Drag
    if 0:
        Cd0 = 0.3

    #######################
    # Ejemplo drag variable
    # CD Coeficiente de drag funncion del Mach

    #Cd=0.125*math.tanh(5*(mach-0.9))+0.425
    #######################
    # Ejemplo lift variable
    # CD Coeficiente de drag funncion del Mach
    #CL_alfa = 0.4*math.tanh(5*(mach-0.9))+0.425

    #######################
    # Ejemplo .308" 168 grain Apendix A pag 217
    if 0:
        # Cd0
        M = [0, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]
        Cd0_exp = [0.14, 0.14, 0.142, 0.16, 0.24, 0.43, 0.449, 0.447, 0.434, 0.41, 0.385, 0.365, 0.35, 0.339, 0.32]

        Cd0 = np.interp(mach, M, Cd0_exp)

        # Cdd2
        Md2 = [0, 0.95, 1, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]
        Cdd2_exp = [2.9, 2.9, 3.0, 3.1, 3.6, 6.5, 7.6, 7.3, 6.8, 6.1, 5.4, 4.4]

        Cdd2 = np.interp(mach, Md2, Cdd2_exp)

        alfa_total2 = ((math.sin(beta))**2 + (math.cos(beta))**2*(math.sin(alfa))**2)
        Cd = Cd0 + Cdd2*alfa_total2

        # Cl_alfa
        Ma = [0, 0.5,  0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5]
        CL_alfa_exp = [1.75, 1.63, 1.45, 1.40, 1.35, 1.3, 1.35, 1.55, 1.7, 1.9, 2.15, 2.32, 2.45, 2.58, 2.68, 2.85]

        CL_alfa = np.interp(mach, Ma, CL_alfa_exp)
        #CL_alfa = CL_alfa * 0.001 # factor de relajacion cl
        CL_alfa = 0.0 # Case B_03, B_04, B_05

        # Cn_p_alfa -> Magnus force
        Cn_p_alfa = -0*0.67

        # Cn_q + Cn_alfa -> Pitching damping debido a q + dot(alfa)
        Cn_q_alfa = -0*0.003

    if globals.Forces_coef_from_txt:
        if not globals.Forces_coef_readed:
            # data = np.loadtxt('Resu_ref/body_coef_txt/bc_baranwonski.txt', delimiter=',', skiprows=3)
            data_raw = np.loadtxt('Resu_ref/body_coef_txt/bc_egipcio.txt', delimiter=',', skiprows=3)
            # data = np.loadtxt('Resu_ref/body_coef_txt/bc_NWU_104pg.txt', delimiter=',', skiprows=3)
            data = data_raw[0:-1]
            globals.data = data
            globals.data_raw = data_raw
            globals.Forces_coef_readed = True

        M = globals.data[:, 0]

        # Drag
        Cd0_exp = globals.data[:, 1]
        Cd0 = np.interp(mach, M, Cd0_exp)

        Cd2_exp = globals.data[:, 2]
        Cd2 = np.interp(mach, M, Cd2_exp)

        alfa_total2 = ((math.sin(beta)) ** 2 + (math.cos(beta)) ** 2 * (math.sin(alfa)) ** 2)
        Cd = Cd0 + Cd2 * alfa_total2

        # Lift
        CL_alfa_exp = globals.data[:, 3]
        CL_alfa = np.interp(mach, M, CL_alfa_exp)

        # Magnus
        Cn_p_alfa_exp = globals.data[:, 4]
        Cn_p_alfa = np.interp(mach, M, Cn_p_alfa_exp)
        Cn_q_alfa = 0

    return [Cd, CL_alfa, Cn_p_alfa, Cn_q_alfa]  # ac'a los toma separados? a los ultimos 2
