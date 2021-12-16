# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import utiles
import math as math

# Grafico la evolucion de los estados

def plot_data(N, t, x, phi, theta, psi):

    fig_size = (12,4)

    fig, axs = plt.subplots(1,3, figsize=fig_size)
    fig.canvas.set_window_title('Coordenadas en marco inercial')
    fig.suptitle('Coordenadas en marco inercial (z=elevacion)')
    for k in range(3):
        axs[k].plot(t, x[:,k])
        axs[k].grid()
    axs[0].set_title('x')
    axs[1].set_title('y')
    axs[2].set_title('z')

    fig, axs = plt.subplots(1,3, figsize=fig_size)
    fig.canvas.set_window_title('Velocidad en marco body')
    fig.suptitle('Velocidad en marco body')
    for k in range(3):
        axs[k].plot(t, x[:,k+3])
        axs[k].grid()
    axs[0].set_title('u')
    axs[1].set_title('v')
    axs[2].set_title('w')

    velocidad_ned = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        velocidad_ned[k] = quat.rotate_vector(xk[3:6])

    fig, axs = plt.subplots(1,4, figsize=fig_size)
    fig.canvas.set_window_title('Velocidad en marco NED')
    fig.suptitle('Velocidad en marco NED')
    for k in range(3):
        axs[k].plot(t, velocidad_ned[:,k])
        axs[k].grid()
    axs[3].plot(t, np.linalg.norm(velocidad_ned,axis=1))
    axs[3].grid()
    axs[0].set_title('vx')
    axs[1].set_title('vy')
    axs[2].set_title('vz')
    axs[3].set_title('mag(V)')

    if 0:
        fig, axs = plt.subplots(1, 4, figsize=fig_size)
        fig.canvas.set_window_title(u'Cuaternión de orientación')
        fig.suptitle(u'Cuaternión de orientación')
        for k in range(4):
            axs[k].plot(t, x[:,k+6])
            axs[k].grid()
        axs[0].set_title('qe')
        axs[1].set_title('qv1')
        axs[2].set_title('qv2')
        axs[3].set_title('qv3')

    if 1:
        fig, axs = plt.subplots(1, 3, figsize=fig_size)
        fig.canvas.set_window_title('Velocidad angular en marco body')
        fig.suptitle('Velocidad angular en marco body')
        for k in range(3):
            axs[k].plot(t, x[:,k+10])
            axs[k].grid()
        axs[0].set_title('p')
        axs[1].set_title('q')
        axs[2].set_title('r')
   # Grafico la evolución de la orientación expresándola mediante los ángulos de Euler


    euler_angles = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        euler_angles[k] = quat.get_euler_anles()

    if 1:
        fig, axs = plt.subplots(1, 3, figsize=fig_size)
        fig.canvas.set_window_title(u'Orientación en ángulos de Euler')
        fig.suptitle(u'Orientación en ángulos de Euler (Deg)')
        for k in range(3):
            axs[k].plot(t, euler_angles[:,k]*180/math.pi)
            axs[k].grid()
        axs[0].set_title('roll (phi)')
        axs[1].set_title('pitch (theta)')
        axs[2].set_title('yaw (psi)')

    if 0:
        #Read results from McCox Fig. 9.2
        xy = np.array([[0],[0]])
        f = open("./Test_Cases/Pitch_yaw_168.dat", "r")
        for line in f:
            frags = line.split()
            xy=np.append(xy, [[float(frags[0])],[float(frags[1])]],axis=1)
        fig, axs = plt.subplots(1, 1, figsize=fig_size)
        fig.canvas.set_window_title(u'pitch vs yaw')
        fig.suptitle(u'2*Clp y 2*Cm_q  y -1*Cm_p_alfa g pitch vs yaw')
        axs.plot(np.rad2deg(euler_angles[:,2]), np.rad2deg(euler_angles[:,1]))
        axs.plot(xy[0,:],xy[1,:], color='green', marker='o', linestyle='dashed',linewidth=1, markersize=2)
        axs.set_title('pitch vs yaw')
        axs.set_xlim(-3, 3)
        axs.set_ylim(-2, 2)
        axs.grid(True)

    if 1:
        fig, axs = plt.subplots(1, 3, figsize=fig_size)
        fig.canvas.set_window_title('Angulos ataque, deslizamiento y total')
        fig.suptitle('Angulos')
        vt = np.linalg.norm(velocidad_ned, axis=1)
        # ángulo de ataque
        alfa = np.arctan2(x[:, 5], x[:, 3])
        alfad = alfa * 180 / math.pi
        axs[0].plot(t, alfad)
        axs[0].grid()
        # ángulo de deslizamiento
        beta = np.arcsin(x[:, 4]/vt)
        betad = beta * 180 / math.pi
        axs[1].plot(t, betad)
        axs[1].grid()
        sin_alfa_t = np.sqrt(((np.sin(beta)) ** 2 + (np.cos(beta)) ** 2 * (np.sin(alfa)) ** 2))
        alfa_t = np.arcsin(sin_alfa_t)
        alfa_t_d = alfa_t * 180 / math.pi
        axs[2].plot(t, alfa_t_d)
        axs[2].grid()
        axs[0].set_title('Alpha')
        axs[1].set_title('Beta')
        axs[2].set_title('Alpha total')

    plt.show(block=False)

    #plt.show()

    return

def plot_force():
    fig_size = (12,4)
    xy = np.array([[0],[0],[0],[0]])
    f = open("./Resultados/Forces_proc.txt", "r")
    line1 = f.readline()
    for line in f:
        frags = line.split(',')
        xy=np.append(xy, [[float(frags[0])],[float(frags[13])],[float(frags[14])],[float(frags[15])]],axis=1)
    xy=np.delete(xy, 0, axis=1)
    fig, axs = plt.subplots(1,3,figsize=fig_size)
    fig.canvas.set_window_title('Forces')

    axs[0].plot(xy[0],xy[1])
    axs[1].plot(xy[0],xy[2])
    axs[2].plot(xy[0],xy[3])

    axs[0].set_title('Force X ')
    axs[1].set_title('Force Y ')
    axs[2].set_title('Force Z ')

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)

    plt.show(block=False)
    #plt.show()

    return


def plot_moments():
    fig_size = (12, 4)
    xy = np.array([[0], [0], [0], [0]])
    f = open("./Resultados/Moments_proc.txt", "r")
    line1 = f.readline()
    for line in f:
        frags = line.split(',')
        xy = np.append(xy, [[float(frags[0])], [float(frags[6])], [float(frags[7])], [float(frags[8])]], axis=1)
    xy = np.delete(xy, 0, axis=1)
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title('Moments')

    axs[0].plot(xy[0], xy[1])
    axs[1].plot(xy[0], xy[2])
    axs[2].plot(xy[0], xy[3])

    axs[0].set_title('Moment X ')
    axs[1].set_title('Moment Y ')
    axs[2].set_title('Moment Z ')

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)

    #plt.show(block=False)
    plt.show()

    return