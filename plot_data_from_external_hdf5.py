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
resu1dir ='./Resultados/ResulRef/'
resu2dir ='./Resultados/Resul_neg_Cm_p_alpha/'
resul = [resu1dir,resu2dir]
leg = ['caso1', 'caso2']
x = []
w = []
vb = []
euler_ang = []
vned = []
t = []
for k in range(np.size(resul)):
    h5f = h5py.File(resul[k]+'Data.hdf5','r')
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
fig_size = (12,4)
fig, axs = plt.subplots(1,3, figsize=fig_size)
fig.canvas.set_window_title('Coordenadas en marco inercial')
fig.suptitle('Coordenadas en marco inercial')
for k in range(s):
    for i in range(r):
        axs[k].plot(t[i], x[i][:,k], label=leg[i])
    axs[k].legend()
    axs[k].grid()
axs[0].set_title('x')
axs[1].set_title('y')
axs[2].set_title('z')

if 0:
    fig3 = plt.figure()
    axs = plt.axes(projection="3d")

    axs.plot3D(x[:,0],x[:,1],x[:,2],'red')
    fig3.canvas.set_window_title('Trayectoria 3D')
    axs.set_title('Trayectoria 3D')
    axs.set_xlabel('x')
    axs.set_ylabel('y')
    axs.set_zlabel('z')

fig, axs = plt.subplots(1,3, figsize=fig_size)
fig.canvas.set_window_title('Velocidad en marco body')
fig.suptitle('Velocidad en marco body')
for k in range(s):
    for i in range(r):
        axs[k].plot(t[i], vb[i][:,k], label=leg[i])
    axs[k].legend()
    axs[k].grid()
axs[0].set_title('u')
axs[1].set_title('v')
axs[2].set_title('w')

fig, axs = plt.subplots(1,4, figsize=fig_size)
fig.canvas.set_window_title('Velocidad en marco NED')
fig.suptitle('Velocidad en marco NED')
for k in range(s):
    for i in range(r):
        axs[k].plot(t[i], vned[i][:,k], label=leg[i])
        axs[3].plot(t[i], np.linalg.norm(vned[i], axis=1))
    axs[k].legend()
    axs[k].grid()
axs[3].grid()
axs[0].set_title('vx')
axs[1].set_title('vy')
axs[2].set_title('vz')
axs[3].set_title('mag(V)')
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
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title('Velocidad angular en marco body')
    fig.suptitle('Velocidad angular en marco body')
    for k in range(s):
        for i in range(r):
            axs[k].plot(t[i], w[i][:,k], label=leg[i])
        axs[k].legend()
        axs[k].grid()
    axs[0].set_title('p')
    axs[1].set_title('q')
    axs[2].set_title('r')

# Grafico la evoluci0n de la orientaci0n expres0ndola mediante los Angulos de Euler
if 1:
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title(u'Orientaci0n en Angulos de Euler')
    fig.suptitle(u'Orientaci0n en 0ngulos de Euler')
    for k in range(s):
        for i in range(r):
            axs[k].plot(t[i], euler_ang[i][:,k]*180/math.pi, label=leg[i])
        axs[k].legend()
        axs[k].grid()
    axs[0].set_title('roll')
    axs[1].set_title('pitch')
    axs[2].set_title('yaw')

    if 1:
        fig, axs = plt.subplots(1, 3, figsize=fig_size)
        fig.canvas.set_window_title('Angulos ataque, deslizamiento y total')
        fig.suptitle('Angulos')
        for i in range(r):
            vt = np.linalg.norm(vned[i], axis=1)
            # angulo de ataque
            alfa = np.arctan2(vb[i][:, 2], vb[i][:, 0])
            alfad = alfa * 180 / math.pi
            axs[0].plot(t[i], alfad)
            # angulo de deslizamiento
            beta = np.arcsin(vb[i][:, 1] / vt)
            betad = beta * 180 / math.pi
            axs[1].plot(t[i], betad)
            sin_alfa_t = np.sqrt(((np.sin(beta)) ** 2 + (np.cos(beta)) ** 2 * (np.sin(alfa)) ** 2))
            alfa_t = np.arcsin(sin_alfa_t)
            alfa_t_d = alfa_t * 180 / math.pi
            axs[2].plot(t[i], alfa_t_d)

        axs[0].grid()
        axs[1].grid()
        axs[2].grid()
        axs[0].set_title('Alpha')
        axs[1].set_title('Beta')
        axs[2].set_title('Alpha total')

def plot_force(resul, leg):
    fig_size = (12,4)
    r = len(resul)
    xy=[]
    for k in range(r):
        tmp = np.array([[0], [0], [0], [0]])
        f = open(resul[k] + 'Forces_proc.txt', 'r')
        line1 = f.readline()
        for line in f:
            frags = line.split(',')
            tmp=np.append(tmp, [[float(frags[0])],[float(frags[13])],[float(frags[14])],[float(frags[15])]],axis=1)
        tmp=np.delete(tmp, 0, axis=1)
        tmp=tmp.transpose()
        xy.append(tmp)
    fig, axs = plt.subplots(1,3,figsize=fig_size)
    fig.canvas.set_window_title('Forces')
    [r, steps, s] = np.shape(xy)

    for i in range(r):
            axs[0].plot(xy[i][:,0],xy[i][:,1], label=leg[i])
            axs[1].plot(xy[i][:,0],xy[i][:,2], label=leg[i])
            axs[2].plot(xy[i][:,0],xy[i][:,3], label=leg[i])

    axs[0].set_title('Force X ')
    axs[1].set_title('Force Y ')
    axs[2].set_title('Force Z ')

    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)

    plt.show(block=False)
    #plt.show()

    return

def plot_moments(resul, leg):
    fig_size = (12, 4)
    r = len(resul)
    xy = []
    for k in range(r):
        tmp = np.array([[0], [0], [0], [0]])
        f = open(resul[k] + "Moments_proc.txt", 'r')
        line1 = f.readline()
        for line in f:
            frags = line.split(',')
            tmp = np.append(tmp, [[float(frags[0])], [float(frags[6])], [float(frags[7])], [float(frags[8])]], axis=1)
        tmp = np.delete(tmp, 0, axis=1)
        tmp = tmp.transpose()
        xy.append(tmp)
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title('Moments')

    for i in range(r):
        axs[0].plot(xy[i][:,0],xy[i][:,1], label=leg[i])
        axs[1].plot(xy[i][:,0],xy[i][:,2], label=leg[i])
        axs[2].plot(xy[i][:,0],xy[i][:,3], label=leg[i])

    axs[0].set_title('Moment X ')
    axs[1].set_title('Moment Y ')
    axs[2].set_title('Moment Z ')

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

