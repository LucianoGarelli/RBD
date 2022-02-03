# script to plot data (Data.hdf5) from external file
# Force_proce must be able in folder "./Resutlados/Forces_proc.txt"
# Moments_proc.txt must be able in folder "./Resutlados/Moments_proc.txt"

import numpy as np
import h5py
import matplotlib.pyplot as plt
import math as math
from mpl_toolkits import mplot3d

## importo hdf5
#h5f = h5py.File('./Resultados/Data.hdf5','r')
resu1dir ='./Resu_ref/Wernert_AIAA2010_7460/Caso_F01/'
resu2dir ='./Resu_ref/Mostafa ASAT-13-FM-03/Caso_E01/'
resul = [resu1dir,resu2dir]
leg = ['caso1', 'caso2']

#plot windows
#ned
fig_size = (12,4)
fig_ned, axs_ned = plt.subplots(1,3, figsize=fig_size)
fig_ned.canvas.set_window_title('Coordenadas en marco inercial')
fig_ned.suptitle('Coordenadas en marco inercial')
#V body
fig_b, axs_b = plt.subplots(1,3, figsize=fig_size)
fig_b.canvas.set_window_title('Velocidad en marco body')
fig_b.suptitle('Velocidad en marco body')
#V ned
fig_n, axs_n = plt.subplots(1,4, figsize=fig_size)
fig_n.canvas.set_window_title('Velocidad en marco NED')
fig_n.suptitle('Velocidad en marco NED')
#Ang veloc
fig_a, axs_a = plt.subplots(1, 3, figsize=fig_size)
fig_a.canvas.set_window_title('Velocidad angular en marco body')
fig_a.suptitle('Velocidad angular en marco body')
#Euler ang
fig_e, axs_e = plt.subplots(1, 3, figsize=fig_size)
fig_e.canvas.set_window_title(u'Orientaci0n en Angulos de Euler')
fig_e.suptitle(u'Orientaci0n en 0ngulos de Euler')
#Angles alpha beta
fig_ang, axs_ang = plt.subplots(1, 3, figsize=fig_size)
fig_ang.canvas.set_window_title('Angulos ataque, deslizamiento y total')
fig_ang.suptitle('Angulos')

for i in range(np.size(resul)):
    x = []
    w = []
    vb = []
    euler_ang = []
    vned = []
    t = []

    h5f = h5py.File(resul[i]+'Data.hdf5','r')
    # parametros hdf5
    x.append(h5f['/Inertial_coord'][:])
    w.append(h5f['/Body_ang_vel'][:])
    vb.append(h5f['/Body_vel'][:])
    euler_ang.append(h5f['/Euler_ang'][:])
    vned.append(h5f['/Inertial_vel'][:])
    t.append(h5f['/Time'][:])
    h5f.close()

# Ploteos
    [r,steps,s] = np.shape(x)

    for k in range(s):
        axs_ned[k].plot(t[0], x[0][:,k], label=leg[i])
        axs_ned[k].legend()
        axs_ned[k].grid()
    axs_ned[0].set_title('x')
    axs_ned[1].set_title('y')
    axs_ned[2].set_title('z')

    if 0:
        fig3 = plt.figure()
        axs = plt.axes(projection="3d")
        axs.plot3D(x[:,0],x[:,1],x[:,2],'red')
        fig3.canvas.set_window_title('Trayectoria 3D')
        axs.set_title('Trayectoria 3D')
        axs.set_xlabel('x')
        axs.set_ylabel('y')
        axs.set_zlabel('z')

    for k in range(s):
        axs_b[k].plot(t[0], vb[0][:,k], label=leg[i])
        axs_b[k].legend()
        axs_b[k].grid()
    axs_b[0].set_title('u')
    axs_b[1].set_title('v')
    axs_b[2].set_title('w')

    for k in range(s):
        axs_n[k].plot(t[0], vned[0][:,k], label=leg[i])
        axs_n[k].legend()
        axs_n[k].grid()
    axs_n[3].plot(t[0], np.linalg.norm(vned[0], axis=1))
    axs_n[3].grid()
    axs_n[0].set_title('vx')
    axs_n[1].set_title('vy')
    axs_n[2].set_title('vz')
    axs_n[3].set_title('mag(V)')
    '''
    # al guardar el Data.hdf5 perdemos el campo de cuaterniones
    if 0:
        fig, axs = plt.subplots(1, 4, figsize=fig_size)
        fig.canvas.set_window_title(u'Cuaterni0n de orientaci0n')
        fig.suptitle(u'Cuaterni0n de orientaci0n')
        for k in range(4):
            axs[k].plot(t, x[:,k+6])
            axs[k].grid()
        axs[0].set_title('qe')
        axs[1].set_title('qv1')
        axs[2].set_title('qv2')
        axs[3].set_title('qv3')
    '''
    if 1:
        for k in range(s):
            axs_a[k].plot(t[0], w[0][:,k], label=leg[i])
            axs_a[k].legend()
            axs_a[k].grid()
        axs_a[0].set_title('p')
        axs_a[1].set_title('q')
        axs_a[2].set_title('r')

# Grafico la evoluci0n de la orientaci0n expres0ndola mediante los Angulos de Euler
    if 1:
        for k in range(s):
            axs_e[k].plot(t[0], euler_ang[0][:,k]*180/math.pi, label=leg[i])
            axs_e[k].legend()
            axs_e[k].grid()
        axs_e[0].set_title('roll')
        axs_e[1].set_title('pitch')
        axs_e[2].set_title('yaw')

    if 1:
        for i in range(r):
            vt = np.linalg.norm(vned[0], axis=1)
            # angulo de ataque
            alfa = np.arctan2(vb[0][:, 2], vb[0][:, 0])
            alfad = alfa * 180 / math.pi
            axs_ang[0].plot(t[0], alfad, label=leg[i])
            # angulo de deslizamiento
            beta = np.arcsin(vb[0][:, 1] / vt)
            betad = beta * 180 / math.pi
            axs_ang[1].plot(t[0], betad, label=leg[i])
            sin_alfa_t = np.sqrt(((np.sin(beta)) ** 2 + (np.cos(beta)) ** 2 * (np.sin(alfa)) ** 2))
            alfa_t = np.arcsin(sin_alfa_t)
            alfa_t_d = alfa_t * 180 / math.pi
            axs_ang[2].plot(t[0], alfa_t_d, label=leg[i])

        axs_ang[0].grid()
        axs_ang[1].grid()
        axs_ang[2].grid()
        axs_ang[0].set_title('Alpha')
        axs_ang[1].set_title('Beta')
        axs_ang[2].set_title('Alpha total')
        axs_ang[0].legend()
        axs_ang[1].legend()
        axs_ang[2].legend()

def plot_force(resul, leg):
    fig_size = (12,4)
    r = len(resul)
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title('Forces')
    for k in range(r):
        xy = []
        tmp = np.array([[0], [0], [0], [0]])
        f = open(resul[k] + 'Forces_proc.txt', 'r')
        line1 = f.readline()
        for line in f:
            frags = line.split(',')
            tmp=np.append(tmp, [[float(frags[0])],[float(frags[13])],[float(frags[14])],[float(frags[15])]],axis=1)
        tmp=np.delete(tmp, 0, axis=1)
        tmp=tmp.transpose()
        xy.append(tmp)
        axs[0].plot(xy[0][:,0],xy[0][:,1], label=leg[k])
        axs[1].plot(xy[0][:,0],xy[0][:,2], label=leg[k])
        axs[2].plot(xy[0][:,0],xy[0][:,3], label=leg[k])


    axs[0].set_title('Force X ')
    axs[0].legend()
    axs[1].set_title('Force Y ')
    axs[1].legend()
    axs[2].set_title('Force Z ')
    axs[2].legend()

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)

    plt.show(block=False)
    #plt.show()

    return

def plot_moments(resul, leg):
    fig_size = (12, 4)
    r = len(resul)
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title('Moments')
    for k in range(r):
        xy = []
        tmp = np.array([[0], [0], [0], [0]])
        f = open(resul[k] + "Moments_proc.txt", 'r')
        line1 = f.readline()
        for line in f:
            frags = line.split(',')
            tmp = np.append(tmp, [[float(frags[0])], [float(frags[6])], [float(frags[7])], [float(frags[8])]], axis=1)
        tmp = np.delete(tmp, 0, axis=1)
        tmp = tmp.transpose()
        xy.append(tmp)
        axs[0].plot(xy[0][:,0],xy[0][:,1], label=leg[k])
        axs[1].plot(xy[0][:,0],xy[0][:,2], label=leg[k])
        axs[2].plot(xy[0][:,0],xy[0][:,3], label=leg[k])

    axs[0].set_title('Moment X ')
    axs[0].legend()
    axs[1].set_title('Moment Y ')
    axs[1].legend()
    axs[2].set_title('Moment Z ')
    axs[2].legend()

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)

    plt.show()

    return

plot_force(resul, leg)
plot_moments(resul, leg)



plt.show()

#print(np.shape(x))
#plt.show()

