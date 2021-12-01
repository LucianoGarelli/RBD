# script to plot data (Data.hdf5) from external file
# Force_proce must be able in folder "./Resutlados/Forces_proc.txt"
# Moments_proc.txt must be able in folder "./Resutlados/Moments_proc.txt"

import numpy as np
import h5py
import matplotlib.pyplot as plt
import math as math
from mpl_toolkits import mplot3d

## importo hdf5
h5f = h5py.File('./Resultados/Data.hdf5','r')
# parametros hdf5
x = h5f['/Inertial_coord'][:] # me importa un vector x de size-> print(np.shape(x)) [1001,3]
w = h5f['/Body_ang_vel'][:]
vb = h5f['/Body_vel'][:]
euler_ang = h5f['/Euler_ang'][:]
vned = h5f['/Inertial_vel'][:]
t = h5f['/Time'][:]
h5f.close()


# Ploteos

fig_size = (12,4)
fig, axs = plt.subplots(1,3, figsize=fig_size)
fig.canvas.set_window_title('Coordenadas en marco inercial')
fig.suptitle('Coordenadas en marco inercial')
for k in range(3):
#    axs[k].plot(t, x[:,k])#,label='Trayectoria en Vacio')    
    axs[k].plot(t, x[:,k])#,label='Trayectoria en Vacio')    
    axs[k].grid()
axs[0].set_title('x')
axs[1].set_title('y')
axs[2].set_title('z')

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
for k in range(3):
#    axs[k].plot(t, x[:,k+3])
    axs[k].plot(t, vb[:,k])
    axs[k].grid()
axs[0].set_title('u')
axs[1].set_title('v')
axs[2].set_title('w')

fig, axs = plt.subplots(1,4, figsize=fig_size)
fig.canvas.set_window_title('Velocidad en marco NED')
fig.suptitle('Velocidad en marco NED')
for k in range(3):
    axs[k].plot(t, vned[:,k])
    axs[k].grid()
axs[3].plot(t, np.linalg.norm(vned,axis=1))
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
    for k in range(3):
        axs[k].plot(t, w[:,k])
        axs[k].grid()
    axs[0].set_title('p')
    axs[1].set_title('q')
    axs[2].set_title('r')

# Grafico la evoluci0n de la orientaci0n expres0ndola mediante los Angulos de Euler
if 1:
    fig, axs = plt.subplots(1, 3, figsize=fig_size)
    fig.canvas.set_window_title(u'Orientaci0n en Angulos de Euler')
    fig.suptitle(u'Orientaci0n en 0ngulos de Euler')
    for k in range(3):
        axs[k].plot(t, euler_ang[:,k]*180/math.pi)
        axs[k].grid()
    axs[0].set_title('roll')
    axs[1].set_title('pitch')
    axs[2].set_title('yaw')

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

    plt.show()

    return

plot_force()
plot_moments()



plt.show()

#print(np.shape(x))
#plt.show()

