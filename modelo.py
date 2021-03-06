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
from moment_coef import moment_coef
from fluid_prop import fluid_prop


m, diam, xcg, ycg, zcg, Ixx, Iyy, Izz, steps, dt = parameters('./Data/data.dat')
# Propiedades fluido
rho, mu , c = fluid_prop(0, 0)
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

    # --- Ecuaciones dinámicas ---#

    # Indica si las fuerzas y momentos se calculan en marco viento o en marco cuerpo
    fuerzas_y_momentos_calculadas_en_marco_body = True

    if fuerzas_y_momentos_calculadas_en_marco_body:
        # por si consideramos viento no nulo. Son las componentes del vector viento en el marco fijo.
        velWind_ned = np.zeros(3)
        # vector vel viento en marco body
        velWind_body = Q_body2ned.T.dot(velWind_ned)
        # velocidad relativa en marco body
        vel_rel_body = x[3:6] - velWind_body
        # supongo que si la velocidad relativa hacia adelante es negativa o nula,
        # no tiene sentido hablar de ángulo de ataque ni de marco viento
        #assert vel_rel_body[0] > 0, "u_vel < 0. %f" % vel_rel_body[0]
        vt = np.linalg.norm(vel_rel_body)
        # ángulo de ataque
        alfa = np.arctan2(vel_rel_body[2], vel_rel_body[0])
        # ángulo de deslizamiento
        beta = np.arcsin(vel_rel_body[1] / vt)

        ca = np.cos(alfa)
        cb = np.cos(beta)
        sa = np.sin(alfa)
        sb = np.sin(beta)
        # matriz de cambio de coordenadas, de marco wind a marco body
        W2B = np.array([[ca * cb, -ca * sb, -sa],
                        [sb, cb, 0],
                        [sa * cb, -sa * sb, ca]])

        # Calculo de coeficientes aerodinamicos eje viento.
        mach = vt/c
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

        # Calculo de aceleraciones en ejes cuerpo
        qdy = 0.5 * rho * vt ** 2

        ####################################
        #Fuerzas en ejes cuerpo

        C_body = np.zeros(6)
        
        C_body[0] = -Cd * math.cos(alfa) * math.cos(beta) + CL_alfa * (1 - math.cos(alfa) ** 2 * math.cos(beta) ** 2)

        C_body[0] = qdy*S*C_body[0]

        C_body[1] = -Cd * math.sin(beta) - CL_alfa * math.cos(alfa) * math.cos(beta) * math.sin(beta)
        # Magnus
        C_body[1] = C_body[1] + (x[10] * diam / (2 * vt)) * Cn_p_alfa * vel_rel_body[2] / vt

        C_body[1] = qdy*S*C_body[1]

        C_body[2] = -Cd * math.sin(alfa) * math.cos(beta) + CL_alfa * math.cos(alfa) * math.sin(alfa) * math.cos(
            beta) ** 2
        # Magnus
        C_body[2] = C_body[2] + (x[10] * diam / (2 * vt)) * Cn_p_alfa * vel_rel_body[1] / vt

        C_body[2] = qdy*S*C_body[2]

        #####################################
        #Momentos ejes cuerpo
        C_body[3] = (x[10] * diam / (2 * vt)) * Clp

        C_body[3] = qdy*S*diam*C_body[3]

        # qt = sqrt(q^2 + r^2) McCoy pag.38 VER Cm_q y Cn_r
        C_body[4] = Cm_alfa * math.sin(alfa) * math.cos(beta) +(x[11] * diam / (2 * vt)) * Cm_q + \
                    (x[10] * diam / (2 * vt)) * Cm_p_alfa * vel_rel_body[1] / vt
        C_body[4] = qdy*S*diam*C_body[4]

        C_body[5] = Cn_beta * math.sin(beta) + (x[12] * diam / (2 * vt)) * Cn_r + \
                    (x[10] * diam / (2 * vt)) * Cm_p_alfa * vel_rel_body[2] / vt

        C_body[5] = qdy*S*diam*C_body[5]

        # antes estaba abajo pero g_body lo pas'e arriba porque necesito para el Forces.txt

        g_body = Q_body2ned.T.dot([0,0,g]) #Q_body2ned[2,:] * g  # multiplico Q_body2ned.T (o sea, la inversa de Q_body2ned) por [0,0,g]^T

        #aca escribir idem arriba pero con Coef[0], Coef[1],....., Coef[5]

        ff = open('./Resultados/Forces.txt', 'ab')
        f_force = np.asarray([dt*(k+1), alfa, beta, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], C_body[0],C_body[1],C_body[2]])
        np.savetxt(ff, [f_force], delimiter=", ", fmt='%1.3e')
        ff.close()

        fm = open('./Resultados/Moments.txt', 'ab')
        m_moment = np.asarray([dt*(k+1), alfa, beta, x[3], x[4], x[5],  C_body[3],C_body[4],C_body[5]])
        np.savetxt(fm, [m_moment], delimiter=", ", fmt='%1.3e')
        fm.close()

        F_body = C_body[0:3]
        M_body = C_body[3:6]

        x_prima[3] = x[12] * x[4] - x[11] * x[5] + g_body[0] + F_body[0] / m
        x_prima[4] = x[10] * x[5] - x[12] * x[3] + g_body[1] + F_body[1] / m
        x_prima[5] = x[11] * x[3] - x[10] * x[4] + g_body[2] + F_body[2] / m

    # Lo siguiente es equivalente a hacer A\b en matlab en vez inv(A)*b (siempre conviene evitar calcular inversas)
    x_prima[10:13] = sp.linalg.solve(inertia_tensor, np.cross(inertia_tensor.dot(x[10:13]),x[10:13]) + M_body, sym_pos=True)
    # sym_pos indica que el tensor de inercia es simétrico y definido positivo, lo que acelera los cálculos


    return x_prima

