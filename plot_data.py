# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import utiles

# Grafico la evolucion de los estados

def plot_data(N, t, x, phi, theta, psi):

    fig_size = (12,4)

    fig, axs = plt.subplots(1,3, figsize=fig_size)
    fig.canvas.set_window_title('Coordenadas en marco inercial')
    fig.suptitle('Coordenadas en marco inercial')
    for k in range(3):
        axs[k].plot(t, x[:,k])
    axs[0].set_title('x')
    axs[1].set_title('y')
    axs[2].set_title('z')

    fig, axs = plt.subplots(1,3, figsize=fig_size)
    fig.canvas.set_window_title('Velocidad en marco body')
    fig.suptitle('Velocidad en marco body')
    for k in range(3):
        axs[k].plot(t, x[:,k+3])
    axs[0].set_title('u')
    axs[1].set_title('v')
    axs[2].set_title('w')

    velocidad_ned = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        velocidad_ned[k] = quat.rotate_vector(xk[3:6])

    fig, axs = plt.subplots(1,3, figsize=fig_size)
    fig.canvas.set_window_title('Velocidad en marco NED')
    fig.suptitle('Velocidad en marco NED')
    for k in range(3):
        axs[k].plot(t, velocidad_ned[:,k])
    axs[0].set_title('vx')
    axs[1].set_title('vy')
    axs[2].set_title('vz')

    if 0:
        fig, axs = plt.subplots(1, 4, figsize=fig_size)
        fig.canvas.set_window_title(u'Cuaternión de orientación')
        fig.suptitle(u'Cuaternión de orientación')
        for k in range(4):
            axs[k].plot(t, x[:,k+6])
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
        fig.suptitle(u'Orientación en ángulos de Euler')
        for k in range(3):
            axs[k].plot(t, euler_angles[:,k])
        axs[0].set_title('roll')
        axs[1].set_title('pitch')
        axs[2].set_title('yaw')

    if 1:
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


    plt.show()

    return
