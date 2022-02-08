# -*- coding: utf-8 -*-
import numpy as np
import os
import h5py as h5py
import utiles

# Grafico la evolucion de los estados

def save_data(N, t, x, Ixx, Iyy, Izz, dir="./Resultados/"):

    CHECK_FOLDER = os.path.isdir(dir)
    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(dir)
        print("Directorio creado: ", dir)

    else:
        print(dir, "Directorio existente.")


    if os.path.isfile(dir + "Data.hdf5"):
        os.remove(dir + "Data.hdf5")
        f = h5py.File(dir + "Data.hdf5", "a")
    else:
        f = h5py.File(dir + "Data.hdf5", "a")

    f.create_dataset('Inertial_coord', data=x[:,0:3])

    f.create_dataset('Time', data=t)

    velocidad_ned = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        velocidad_ned[k] = quat.rotate_vector(xk[3:6])
    f.create_dataset('Inertial_vel', data=velocidad_ned)


    f.create_dataset('Body_vel', data=x[:,3:6])
    f.create_dataset('Body_ang_vel', data=x[:,10:13],dtype=np.float64) # modificamos formato de numero | np.float64

    velocidad_ang_ned = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        velocidad_ang_ned[k] = quat.rotate_vector(xk[10:13])
    f.create_dataset('Inertial_ang_vel', data=velocidad_ang_ned,dtype=np.float64) # modificamos formato de numero)

    euler_angles = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        euler_angles[k] = quat.get_euler_anles()
    f.create_dataset('Euler_ang', data=euler_angles)

    ##matriz_rot = np.empty((N+1,3))
    #for k, xk in enumerate(x):
    #    matriz_rot[k] = quat.rotate_vector()
    #f.create_dataset('matriz_rot', data=euler_angles)
    
    f.close()

    return velocidad_ang_ned
