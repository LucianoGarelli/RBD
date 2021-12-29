# -*- coding: utf-8 -*-
"""
@author: lgenzelis
"""
import math as math
import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt
import modelo
import utiles
import conservation as cons
import os

from parameters import parameters
from fluid_prop import fluid_prop
from initial_cond import initial_cond
import plot_data as plt_data
import save_data as sv
from proce import proce


def main():

    print '#########################'
    print 'Reading parameters.'
    np.set_printoptions(precision=3)
    np.set_printoptions(suppress=True)
    # Parametros proyectil
    m, diam, xcg, ycg, zcg, Ixx, Iyy, Izz, steps, dt = parameters('./Data/data.dat')
    print "Masa=", m, "[kg]"
    print "Diam=", diam, "[m]"
    print "xcg, ycg, zcg=", xcg, ycg, zcg, "[m]"
    print "Ixx, Iyy, Izz=", Ixx, Iyy, Izz, "[m^4]"
    print "Steps, Dt=", steps, dt

    print '#########################'
    print 'Reading fluid properties.'
    # Propiedades fluido
    rho, mu, c = fluid_prop(0, 0)
    print "Densidad=", rho, "[kg/m3]"
    print "mu=", mu, "[kg/m s]"
    print "c=", c, "[m/s]"

    print '#########################'
    print 'Reading Initial Conditions.'
    # Condiciones iniciales
    V, alfa, beta, p, q, r, phi, theta, psi, XE, YE, ZE = initial_cond('./Data/Initial_cond.dat')
    print "Velocity=", V, "[m/s]"
    print "alfa, beta=", alfa, beta, "[deg]"
    print "p,q,r=", p, q, r, "[RPM]"
    print "phi, theta, psi=", phi, theta, psi, "[deg]"
    print "XE, YE, ZE=", XE, YE, ZE, "[m]"

    MYDIR = ("Resultados")
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("Directorio Resultados creado: ", MYDIR)

    else:
        print(MYDIR, "Directorio Resultados existente.")

    #File to write force coeff
    ff = open("./Resultados/Force_coef.txt", "w")  # xq le pasamos los las fuerzas de todo el CFD??
    ff.write(" # Time,   Mach,     alfa,     beta,     Cd,     CL_alfa,     Cn_p_alfa,     Cn_q_alfa \n")
    ff.close()
    #File to write moment coeff
    fm = open("./Resultados/Moment_coef.txt", "w")
    fm.write(" # Time,   Mach,     alfa,     beta,     Clp,     Cm_alfa,     Cm_p_alfa,     Cm_q,     Cn_beta,"
             "     Cn_r \n")
    fm.close()

    # ---------------------------------------------
    # File to write forces - Body Frame
    ff = open("./Resultados/Forces.txt", "w")
    ff.write("# Time,     alpha,     beta,     V_inf,    u(v_body_X), v(v_body_Y), w(v_body_Z),  p, "
             "      q,         r,         gx,         gy,       gz,      FX_body,    FY_body,    FZ_body  \n")
    ff.close()
    # File to write moment - Body Frame
    fm = open("./Resultados/Moments.txt", "w")
    fm.write("# Time,     alpha,      beta,     p,      q,      r,    MX,     MY,     MZ \n")
    fm.close()
    #---------------------------------------------
    
    # Convierto de deg->rad y de RPM->rad/s
    alfa = np.deg2rad(alfa)
    beta = np.deg2rad(beta)

    phi = np.deg2rad(phi)
    theta = np.deg2rad(theta)
    psi = np.deg2rad(psi)

    p = p * 2 * math.pi / 60
    q = q * 2 * math.pi / 60
    r = r * 2 * math.pi / 60

    Ts = dt  # periodo de discretización
    N = steps # el tiempo total de la simulación es N*Ts

    # Estado inicial
    x0 = np.zeros(13)
    q0 = utiles.Quaternion.fromEulerAngles(roll=phi, pitch=theta, yaw=psi)  # orientación inicial, obtenida a partir de los ángulos de Euler
    x0[6] = q0.d
    x0[7:10] = q0.v
    x0[0:3] += [XE, YE, ZE]
    x0[3] = V*math.cos(alfa)*math.cos(beta) # velocidad inicial en body
    x0[4] = V*math.sin(beta)
    x0[5] = V*math.sin(alfa)*math.cos(beta)
    x0[10] = p # rotacion inicial en body
    x0[11] = q
    x0[12] = r

    x = np.zeros((N+1,13))
    x[0] = x0
    u = np.zeros((N,4))  # tomo 4 aciones de control, por ejemplo. El control correspondiente a u[k] se mantiene constante desde el instante k*Ts hasta el instante (k+1)*Ts
    t = np.arange(0, N+1) * Ts
    for k in range(N):
        x[k + 1] = sp.integrate.odeint(lambda _x, _t: modelo.ED_cuaterniones(_x, u[k], k, _t), x[k], [k*Ts, (k+1)*Ts],
                                       rtol=1e-6, atol=1e-6)[-1]

        #output = sp.integrate.odeint(lambda _x, _t: modelo.ED_cuaterniones(_x, u[k], k, _t), x[k], [k * Ts, (k + 1) * Ts],
        #                               rtol=1e-12, atol=1e-12, full_output=True)
        #x[k+1] = output[0][-1]
        #print np.unique(output[1]["tsw"])
        #x[k + 1] = solucion[-1]

        #Renormalizacion quaterniones
        #x[k+1,6:10] /= np.linalg.norm(x[k+1,6:10])
        t_N = k+1
        if x[1,2] < -1.0:
        #15 yardas
        #if solucion[1,0] >= 13.7:
            break
        if (k % 5000) == 0:
            sv.save_data(N,t,x,Ixx,Iyy,Izz)
            proce('./Resultados/Forces.txt')
            proce('./Resultados/Moments.txt')


    print '#########################'
    print 'Resultados Generales'
    print "Altura maxima =", np.amax(x[:,2])  , "[m]"
    print "Distancia =", math.sqrt((x[t_N,0]**2+x[t_N,1]**2)) , "[m]"
    print "Altura para tiempo final =", t[t_N], "[s]", (x[t_N,2]) , "[m]"

    quat = utiles.Quaternion(x[t_N,6], x[t_N,7:10])
    velocidad_ned = quat.rotate_vector(x[t_N,3:6])

    print "Velocidades para tiempo final Vx,Vy,Vz=", velocidad_ned , "[m/s]"
    print "Coordenadas para tiempo final Xned,Yned,Zned=", x[t_N,0:3] , "[m]"

    sv.save_data(N,t,x,Ixx,Iyy,Izz)
    #cons.conservation(N, x, Ixx, Iyy, Izz)
    proce('./Resultados/Force_coef.txt')
    proce('./Resultados/Forces.txt')
    proce('./Resultados/Moment_coef.txt')
    proce('./Resultados/Moments.txt')
    plt_data.plot_data(N, t, x, phi, theta, psi)
    plt_data.plot_force()
    plt_data.plot_moments()


if __name__ == "__main__":
    main()
