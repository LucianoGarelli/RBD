# -*- coding: utf-8 -*-
"""
@author: lgenzelis
"""

import numpy as np
import scipy as sp
import utiles
import math as math

from parameters import parameters
from force_coef import force_coef
from force_coef_body import force_coef_body
from moment_coef import moment_coef
from moment_coef_body import moment_coef_body
from fluid_prop import fluid_prop


m, diam, xcg, ycg, zcg, Ixx, Iyy, Izz, steps, dt = parameters('./Data/data.dat')

g= 9.81  # aceleración de la gravedad

inertia_tensor = np.array([[Ixx, 0., 0.],  # el tensor de inercia en el marco body
                           [0., Iyy, 0.],
                           [0., 0., Izz]])
# Superficie Ref
S = math.pi * (0.5 * diam) ** 2

def ED_cuaterniones(x, u, k, t):
    '''
    :param x: Vector de estados
    :param u: Vector de acciones de control => podria incidir incidir sobre las fins
    :return: Vector de derivadas de los estados
    x[0] -> x_ned (north)
    x[1] -> y_ned (east)
    x[2] -> h_ned = -z_ned (z_ned = down)
    x[3] -> u (vel x en marco body)
    x[4] -> v (vel y en marco body)
    x[5] -> w (vel z en marco body)
    x[6] -> q_e (parte escalar del cuaternión de orientación)
    x[7] -> q_v1 (primera componente de la parte vectorial del cuaternión de orientación)
    x[8] -> q_v2 (segunda componente de la parte vectorial del cuaternión de orientación)
    x[9] -> q_v3 (tercera componente de la parte vectorial del cuaternión de orientación)
    x[10] -> p (1era componente de la velocidad angular en marco body)
    x[11] -> q (2da componente de la velocidad angular en marco body)
    x[12] -> r (3era componente de la velocidad angular en marco body)
    x[13] -> alfa
    x[14] -> beta
    '''
    # t = tiempo absoluto
    x_prima = np.zeros_like(x)

    q_body2ned = utiles.Quaternion(x[6], x[7:10])
    Q_body2ned = q_body2ned.calc_rot_matrix()

    #--- Ecuaciones cinemáticas ---#

    x_prima[0:3] = Q_body2ned.dot(x[3:6])  # paso la velocidad de body[3:6] a ned[0:3]
    x_prima[2] *= -1. # para tener altura en vez de "profundidad"

    q_prima = q_body2ned.mult_cuat_times_vec(x[10:13]*.5)
    x_prima[6] = q_prima.d
    x_prima[7:10] = q_prima.v

    # Propiedades fluido
    rho, mu, c = fluid_prop(x[2], 0)
    # --- Ecuaciones dinámicas ---#

    # Indica si las fuerzas y momentos se calculan en marco viento / marco cuerpo / marco ned
    fuerzas_y_momentos_calculadas_en_marco_body = True

    # por si consideramos viento no nulo. Son las componentes del vector viento en el marco fijo.
    velWind_ned = np.zeros(3)
    # vector vel viento en marco body
    velWind_body = Q_body2ned.T.dot(velWind_ned)
    # velocidad relativa en marco body
    vel_rel_body = x[3:6] - velWind_body
    #vel_rel_ned = x[0:3] - velWind_ned
    # supongo que si la velocidad relativa hacia adelante es negativa o nula,
    # no tiene sentido hablar de ángulo de ataque ni de marco viento
    # assert vel_rel_body[0] > 0, "u_vel < 0. %f" % vel_rel_body[0]
    vt = np.linalg.norm(vel_rel_body)
    # ángulo de ataque
    alfa = np.arctan2(vel_rel_body[2], vel_rel_body[0])
    alfad = alfa * 180 / math.pi
    # ángulo de deslizamiento
    beta = np.arcsin(vel_rel_body[1] / vt)
    betad = beta * 180 / math.pi

    sin_alfa_t = math.sqrt(((math.sin(beta)) ** 2 + (math.cos(beta)) ** 2 * (math.sin(alfa)) ** 2))
    alfa_t = math.asin(math.sqrt(((math.sin(beta)) ** 2 + (math.cos(beta)) ** 2 * (math.sin(alfa)) ** 2)))

    ca = np.cos(alfa)
    cb = np.cos(beta)
    sa = np.sin(alfa)
    sb = np.sin(beta)

    mach = vt / c

    if fuerzas_y_momentos_calculadas_en_marco_body:

        # matriz de cambio de coordenadas, de marco wind a marco body
        W2B = np.array([[ca * cb, -ca * sb, -sa],
                        [sb, cb, 0],
                        [sa * cb, -sa * sb, ca]])

        Cx0, Cx2, Cna, Cypa = force_coef_body(mach,alfa,beta)
        Cma, Cmq, Clp, Cnpa = moment_coef_body(mach,alfa,beta)

        ff = open('./Resultados/Force_coef.txt', 'ab')
        f_coef = np.asarray([dt*(k +1),  mach, alfa, beta, Cx0, Cx2, Cna])
        np.savetxt(ff, [f_coef], delimiter=", ", fmt='%1.3e')
        ff.close()

        fm = open('./Resultados/Moment_coef.txt', 'ab')
        m_coef = np.asarray([dt*(k+1), mach, alfa, beta, Cma, Cmq, Clp])
        np.savetxt(fm, [m_coef], delimiter=", ", fmt='%1.3e')
        fm.close()

        # Calculo presión dinámica
        qdy = 0.5 * rho * vt ** 2

        ####################################
        #Fuerzas en ejes cuerpo
        C_body = np.zeros(6)
        NED_forces = np.zeros(6)
        u = vel_rel_body[0]
        v = vel_rel_body[1]
        w = vel_rel_body[2]
        C_body[0] = -qdy*S*(Cx0 + Cx2*(v**2 + w**2)/vt**2)
        C_body[1] = -qdy*S*(Cna*v/vt - Cypa*(x[10] * diam / (2 * vt))*(w/vt))
        C_body[2] = -qdy*S*(Cna*w/vt + Cypa*(x[10] * diam / (2 * vt))*(v/vt))

        #####################################
        #Momentos ejes cuerpo
        C_body[3] = qdy*S*diam*(x[10] * diam / (2 * vt)) * Clp
        # qt = sqrt(q^2 + r^2) McCoy pag.38 VER Cm_q y Cn_r
        C_body[4] = qdy*S*diam*(Cma*w/vt + Cmq*(x[11] * diam / (2 * vt)) + Cnpa*(x[10] * diam / (2 * vt))*(v/vt))
        C_body[5] = qdy*S*diam*(-Cma*v/vt + Cmq*(x[12] * diam / (2 * vt)) + Cnpa*(x[10] * diam / (2 * vt))*(w/vt))

        #
        # ver g_body y NED_forces esto estaba en C_body
        # antes estaba abajo pero g_body lo pas'e arriba porque necesito para el Forces.txt
        g_body = Q_body2ned.T.dot([0,0,g]) #Q_body2ned[2,:] * g  # multiplico Q_body2ned.T (o sea, la inversa de Q_body2ned) por [0,0,g]^T
        NED_forces[0:3] = Q_body2ned.dot(C_body[0:3])

        ff = open('./Resultados/Forces.txt', 'ab')
        f_force = np.asarray([dt*(k+1), alfad, betad, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], C_body[0],C_body[1],C_body[2]])
        np.savetxt(ff, [f_force], delimiter=", ", fmt='%1.3e')
        ff.close()

        fm = open('./Resultados/Moments.txt', 'ab')
        m_moment = np.asarray([dt*(k+1), alfad, betad, x[3], x[4], x[5], C_body[3], C_body[4], C_body[5]])
        np.savetxt(fm, [m_moment], delimiter=", ", fmt='%1.3e')
        fm.close()

        Forces = C_body[0:3]
        Moments = C_body[3:6]

        #aca escribir idem arriba pero con Coef[0], Coef[1],....., Coef[5]
    else:
        # Coeficientes aerodinamicos eje viento
        Cd, CL_alfa, Cn_p_alfa, Cn_q_alfa = force_coef(mach, alfa, beta)
        Clp, Cm_alfa, Cm_p_alfa, Cm_q, Cn_beta, Cn_r = moment_coef(mach, alfa, beta)
        
        ff = open('./Resultados/Force_coef.txt', 'ab')
        f_coef = np.asarray([dt*(k +1),  mach, alfa, beta, Cd, CL_alfa, Cn_p_alfa, Cn_q_alfa])
        np.savetxt(ff, [f_coef], delimiter=", ", fmt='%1.3e')
        ff.close()

        fm = open('./Resultados/Moment_coef.txt', 'ab')
        m_coef = np.asarray([dt*(k+1), mach, alfa, beta, Clp, Cm_alfa, Cm_p_alfa, Cm_q, Cn_beta, Cn_r])
        np.savetxt(fm, [m_coef], delimiter=", ", fmt='%1.3e')
        fm.close()
            
        # Calculo presión dinámica
        qdy = 0.5 * rho * vt ** 2

        ####################################
        #Fuerzas en ejes viento
        C_Wind = np.zeros(6)
        NED_forces = np.zeros(6)
        u = vel_rel_body[0]
        v = vel_rel_body[1]
        w = vel_rel_body[2]
        C_Wind[0] = -qdy*S*Cd*ca*cb + qdy*S*CL_alfa*(sb**2 + sb**2 * cb**2 ) + 0
        C_Wind[1] = -qdy*S*Cd*sb + qdy*S*CL_alfa*(-ca*sb*cb) + qdy*S*Cn_p_alfa*x[10]*(-sa*cb)
        C_Wind[2] = -qdy*S*Cd*sa*cb + qdy*S*CL_alfa*(-sa*ca*cb**2) + qdy*S*Cn_p_alfa*x[10]*(sb)
            
        #####################################
        #Momentos ejes viento
        C_Wind[3] = 0 + 0 + qdy*S*diam*(diam/vt)*Clp*x[10] + 0 
        # qt = sqrt(q^2 + r^2) McCoy pag.38 VER Cm_q y Cn_r
        C_Wind[4] = qdy*S*diam*Cm_alfa*(sa*cb) + qdy*S*diam*(diam/vt)*Cm_p_alfa*(-sb) + 0 + qdy*S*diam*(diam/vt)*Cm_q * x[11]
        C_Wind[5] = qdy*S*diam*Cm_alfa*(-sb) + qdy*S*diam*(diam/vt)*Cm_p_alfa*(-sa*cb) + 0 + qdy*S*diam*(diam/vt)*Cm_q * x[12]
            
        # ver como plantear la gravedad en ejes viento y eso de NED_forces
        # g_body no lo cambié, la fórmula es similar al g_body
        g_body = Q_body2ned.T.dot([0,0,g]) #Q_body2ned[2,:] * g  # multiplico Q_body2ned.T (o sea, la inversa de Q_body2ned) por [0,0,g]^T
        #NED_forces[0:3] = Q_body2ned.dot(C_body[0:3])

        ff = open('./Resultados/Forces.txt', 'ab')
        f_force = np.asarray([dt*(k+1), alfad, betad, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], C_Wind[0],C_Wind[1],C_Wind[2]])
        #f_force = np.asarray([dt*(k+1), alfad, betad, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], NED_forces[0],NED_forces[1],NED_forces[2]])
        np.savetxt(ff, [f_force], delimiter=", ", fmt='%1.3e')
        ff.close()

        fm = open('./Resultados/Moments.txt', 'ab')
        # modificar este para momentos wind
        m_moment = np.asarray([dt*(k+1), alfad, betad, x[3], x[4], x[5], C_Wind[3], C_Wind[4], C_Wind[5]])
        np.savetxt(fm, [m_moment], delimiter=", ", fmt='%1.3e')
        fm.close()

        Forces = C_Wind[0:3]
        Moments = C_Wind[3:6]

        #
        # ver que hacemos con g_body y NED_forces, antes estaban arriba
        #

    x_prima[3] = x[12] * x[4] - x[11] * x[5] + g_body[0] + Forces[0] / m
    x_prima[4] = x[10] * x[5] - x[12] * x[3] + g_body[1] + Forces[1] / m
    x_prima[5] = x[11] * x[3] - x[10] * x[4] + g_body[2] + Forces[2] / m

    # Lo siguiente es equivalente a hacer A\b en matlab en vez inv(A)*b (siempre conviene evitar calcular inversas)
    x_prima[10:13] = sp.linalg.solve(inertia_tensor, np.cross(inertia_tensor.dot(x[10:13]), x[10:13]) + Moments, sym_pos=True)
    # sym_pos indica que el tensor de inercia es simétrico y definido positivo, lo que acelera los cálculos



    return x_prima

