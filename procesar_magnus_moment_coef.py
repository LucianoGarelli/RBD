# funci'on de prueba para procesar magnus coef (momento)
import numpy as np
import globals
def procesar_magnus_moment_coef():
    # colocar en el vecto n_magnus los distintos 'angulos que se disponen de coeficientes

    # cant de angulos que se evaluan los coef
    n_magnus = [0.0, 2.0, 5.0, 10.0] # egipcio.txt
    #n_magnus = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0] # bc_baranwonski
    #n_magnus = [1] # NWU_104pg
    n = np.size(n_magnus)
    
    #data = np.loadtxt('./bc_NWU_104pg.txt', delimiter=',', skiprows=3)
    #data = np.loadtxt('./delete_bc_baranwonski.txt', delimiter=',', skiprows=3)
    #data = np.loadtxt('./delete_bc_egipcio.txt', delimiter=',', skiprows=3)
    [rows,cols] = np.shape(globals.data)
    tot = rows*n
    M_0 = globals.data[:,0]


    # Mach para doble interpolaci'on
    globals.Mpa = np.zeros(tot)
    for i in range(rows):
        for k in range(np.size(n_magnus)):
            globals.Mpa[i*n+k] = M_0[i]
        print('globals.Mpa')
        print(globals.Mpa)
        print(np.shape(globals.Mpa))


    # angulos de magnus
    globals.ang = []
    for ii in range(rows):
        globals.ang = np.append(globals.ang,n_magnus)
    print('globals.ang')
    print(globals.ang)
    print(np.shape(globals.ang))


    # Cnpa
    globals.Cnpa_proce = np.zeros(tot)
    for i in range(rows):
        for k in range(np.size(n_magnus)):
            globals.Cnpa_proce[i*n+k] = globals.data[i,k+8]
    print('globals.Cnpa_proce')
    print(globals.Cnpa_proce)
    print(np.shape(globals.Cnpa_proce))


    return []
