# funci'on de prueba para procesar magnus coef (momento)
import numpy as np
import globals
def procesar_magnus_moment_coef():
    # colocar en el vecto n_magnus los distintos 'angulos que se disponen de coeficientes

    [rows,cols] = np.shape(globals.data)

    n = cols - 8
    n_magnus = globals.data [rows-1,8:8+n]
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
