# script for post-processing results from MHE estimation and fits data as body/wind_coef for inputs for the RBD
import numpy as np
import matplotlib.pyplot as plt
import os

dir_load = '/home/ntrivi/Documents/Tesis/Estim_MHE/Resultados/'
file_name = 'Coef_estim.txt'
data = np.loadtxt(dir_load + file_name, delimiter=',', skiprows=1)
# characteristic of the estimation to incorporate in the txt file
Case = 'F01'
Nt = '1'
MHE = 'MPC casadi tools' # MPC tools or ad-hod
Noise = 'True'

x = data
#  Time  ,    Mach  ,  alpha_tot,    Cd0   ,  Cl_alpha,    Cd2   , Cn_p_alpha,    Clp    ,  Cm_alpha , Cm_p_alpha,    Cm_q \n

n1 = np.size(x[:,1])
[n1,n2] = np.shape(x)
x0 = 3 # number of data columns written in Coef_estim.txt: time, mach, alpha_tot
ncoef = n2-3 # number of coef obtained in Estim_MHE

# Mach array to organice the coef
M = np.array([0.01, 0.60, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.35, 1.50, 1.75, 2.00, 5.00])
n3 = np.size(M)
aref = np.array([0.0, 2.0, 5.0, 10.0]) # alpha ref
    
acum = np.zeros((n3,ncoef))
cont = np.zeros((n3,ncoef))
coef = np.zeros((n3,ncoef))
i = 0 # inicializamos el contador

for j in range(n1):
    for i in range(n3):
        for k in range(3,ncoef+x0): # vectorizado en funcion de N coef
            if M[i] <= x[j,1] and x[j,1] < M[i+1]:
                if x[j,1] > (M[i+1]+M[i])/2:
                    acum[i+1,k-x0] += x[j,k]
                    cont[i+1,k-x0] += 1
                else:
                    acum[i,k-x0] += x[j,k]
                    cont[i,k-x0] += 1

for k in range(n3):
    for j in range(ncoef):
        if cont[k,j] != 0:
            coef[k,j] = acum[k,j]/cont[k,j]

# Subsonic regimen adopted constant
coef[0,:] = coef[1,:]

# FIXME ! se puede mejorar con un vector de re-ordenamiento
coef_reorg = np.zeros((n3,ncoef))
coef_reorg[:,0] = coef[:,0] # Cd0
coef_reorg[:,1] = coef[:,2] # Cd2
coef_reorg[:,2] = coef[:,1] # Cla
coef_reorg[:,3] = coef[:,3] # Cnypa
coef_reorg[:,4] = coef[:,4] # Clp
coef_reorg[:,5] = coef[:,5] # Cma
coef_reorg[:,6] = coef[:,7] # Cmq
coef_reorg[:,7] = coef[:,6] # Cnpa

M = M.reshape(-1,1)
coef2 = np.concatenate((M,coef_reorg),axis=1)

# File to write force coeff
dir_load_w = '/home/ntrivi/Documents/Tesis/RBD/Resu_ref/body_coef_txt/'
file_name_w = 'Coef_estim_from_MHE_' + Case + '_Nt_' + Nt + '_Noise_' + Noise +'.txt'
file_w = dir_load_w + file_name_w

ff = open(file_w, "w")
ff.write("# Coef readed from Estimation_MHE\n")
ff.write("# re-group data fitting for RBD body_coef from " + dir_load + file_name + ', Case: ' + Case + ', Nt: ' + Nt + ', MHE method: ' + MHE + ', Noise added: ' + Noise + "\n")# 3 lines left because the other script skip the first 3 lines
ff.write("# Mach,       Cd0,       Cd2,    Cl_alpha,  Cn_p_alpha,   Clp,    Cm_alpha,    Cm_q ,   Cm_p_alpha\n")
np.savetxt(ff, coef2, delimiter=", ", fmt='%1.3e')
ff.close()

'''
organization for the double entry-table
#------------------------------------------
aref = np.array([0.0, 2.0, 5.0, 10.0])  # alpha ref
aref = np.radians(aref)
#pp = foo(x[:,2],aref)

cmpa = x[:,-1]

n1 = np.size(cmpa)

# Mach array to organice the coef
M = np.array([0.01, 0.60, 0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.35, 1.50, 1.75, 2.00, 5.00])
n3 = np.size(M)

# for alpha ref
acum_a = np.zeros(n2)
cont_a = np.zeros(n2)
coef_a = np.zeros(n2)
acum_cm = np.zeros(n2)
cont_cm = np.zeros(n2)
coef_cm = np.zeros(n2)
i = 0 # inicializamos el contador

for j in range(n1):
    p = cmpa[j]
    for i in range(n2):
        if aref[i] <= x[j] and x[j] < aref[i+1]:
            if x[j] > (aref[i+1]+aref[i])/2:
                acum_a[i+1] += x2[j]
                cont_a[i+1] += 1
                acum_cm[i+1] += x2[j]
                cont_cm[i+1] += 1
            else:
                acum_a[i] += x2[j]
                cont_a[i] += 1
                acum_cm[i] += x2[j]
                cont_cm[i] += 1

print('Verification',bool(sum(cont)==n1))

for k in range(np.size(acum)):
    if cont_a[k] != 0:
        coef_a[k] = acum[k]/cont[k]
        coef_cm[k] = acum[k]/cont[k]

#termina alpha ref
#------------------------------------------
#'''
print('\n*-----------------------------------------------*\n')
print('FIN! - OK \n')
