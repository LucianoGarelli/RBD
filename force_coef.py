import math as math
import numpy as np

def force_coef(mach, alfa, beta):
    #######################
    # Datos Ejemplo Carlucci pag. 244
    # CD Coeficiente de Drag
    # #Cd0 = 0.3
    #######################

    #######################
    # Ejemplo .308" 168 grain Apendix A pag 217

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


    # Cn_p_alfa -> Magnus force
    Cn_p_alfa = -0*0.67

    # Cn_q + Cn_alfa -> Pitching damping debido a q + dot(alfa)
    Cn_q_alfa = -0*0.003

    return [Cd, CL_alfa, Cn_p_alfa, Cn_q_alfa]  # ac'a los toma separados? a los ultimos 2
