# -*- coding: utf-8 -*-
# funcion que verifica la conservacion de momento sobre eje cuerpo y eje tierra
import numpy as np
import os
import h5py as h5py
import utiles
import math as math
import numpy as np


def conservation(N, x, Ixx, Iyy, Izz):

    print '#########################' # señal de que esta bien hecha la funcion
    print "entro en conservation function"
        
    #   print(np.finfo(float).eps)
    eps = np.finfo(float).eps
    print "eps = ", eps

    L1b = Ixx * x[:,10]
    L2b = Iyy * x[:,11]
    L2b += eps 
    L3b = Izz * x[:,12]
    L3b += eps
    
    Lb = L1b + L2b + L3b;
    Ltotb = L1b ** 2 + L2b ** 2 + L3b ** 2

    velocidad_ang_ned = np.empty((N+1,3))
    for k, xk in enumerate(x):
        quat = utiles.Quaternion(xk[6], xk[7:10])
        velocidad_ang_ned[k] = quat.rotate_vector(xk[10:13])

    L1ned = Ixx * velocidad_ang_ned[:,0]
    L2ned = Iyy * velocidad_ang_ned[:,1]
    L3ned = Izz * velocidad_ang_ned[:,2]

    Ltotned = L1ned ** 2 + L2ned ** 2 + L3ned ** 2

    print "max(L2b)",max(L2b),"min(L2b)",min(L2b)
    print "max(L3b)",max(L3b),"min(L3b)",min(L3b)
    print '#########################'

    if 1:
        print "Calculando la conservacion como - abs(max(L1b) / min(L1b)) - 1:"
        print "Conservacion Momento Cuerpo L1b: =", abs(max(L1b) / min(L1b)) - 1
        print "Conservacion Momento Cuerpo L2b: =", abs(max(L2b) / min(L2b)) - 1
        print "Conservacion Momento Cuerpo L3b: =", abs(max(L3b) / min(L3b)) - 1
        print "Conservacion Momento Cuerpo Ltotb: =", abs(max(Ltotb) / min(Ltotb)) - 1
        print 'Conservacion Momento Tierra: L1ned: =',abs(max(L1ned)/min(L1ned))-1
        print 'Conservacion Momento Tierra: L2ned: =',abs(max(L2ned)/min(L2ned))-1
        print 'Conservacion Momento Tierra: L3ned: =',abs(max(L3ned)/min(L3ned))-1
        print 'Conservacion Momento Tierra: Ltotned: =',abs(max(Ltotned)/min(Ltotned))-1
        print "Ixx =",Ixx," Iyy =",Iyy," Izz =",Izz

    if 0:
        print "Calculando la conservacion como - max(L2b)-min(L2b):"
        print "Conservacion Momento Cuerpo L1b: =", max(L1b)-min(L1b)
        print "Conservacion Momento Cuerpo L2b: =", max(L2b)-min(L2b)
        print "Conservacion Momento Cuerpo L3b: =", max(L3b)-min(L3b)
        print "Conservacion Momento Cuerpo L3tot: =", max(Ltotb)-min(Ltotb)
        print "Conservacion Momento Cuerpo L1ned: =", max(L1ned)-min(L1ned)
        print "Conservacion Momento Cuerpo L2ned: =", max(L2ned)-min(L2ned)
        print "Conservacion Momento Cuerpo L3ned: =", max(L3ned)-min(L3ned)
        print "Conservacion Momento Cuerpo L3totned: =", max(Ltotned)-min(Ltotned)
    

    Inertia = np.array([Ixx, Iyy, Izz])
    min_inertia = np.min(Inertia)
    print "Inertia_min", min_inertia
    Inertia_rel = Inertia/min_inertia

    print "Inertia_rel = ", Inertia_rel

    return
