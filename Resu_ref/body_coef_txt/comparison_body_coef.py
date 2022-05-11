# archivo prueba leer data
import numpy as np
import matplotlib.pyplot as plt
				
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Carga coeficientes seg'un distintas bibliograf'ias
data_eg = np.loadtxt('bc_egipcio.txt', delimiter=',', skiprows=3)
data_eg = data_eg[0:len(data_eg[:,0])-1,:]
# Encabezado del txt:
# M	Ca	Cad2	Cna	Cypa	Clp	Cma	Cmq	Cnpa			
#								0	2	5	10
data_ba = np.loadtxt('bc_baranwonski.txt', delimiter=',', skiprows=3)
data_ba = data_ba[0:len(data_ba[:,0])-1,:]
# Encabezado del txt:
# M	Ca	Cad2	Cna	Cypa	Clp	Cma	Cmq	Cnpa			
#								0	1       2	3      5	10       20
data_NWU = np.loadtxt('bc_NWU_104pg.txt', delimiter=',', skiprows=3)
data_NWU = data_NWU[0:len(data_NWU[:,0])-1,:]
# M	Ca	Cad2	Cna	Cypa	Clp	Cma	Cmq	Cnpa	
# aca el magnus moment tiene un solo valor y no indica el 'angulo									
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# characteristic of the estimation to incorporate in the txt file MHE2RBD
Case = 'F01'
Nt = '1'
MHE = 'MPC casadi tools' # MPC tools or ad-hod
Noise = 'True'

file_name_data_MHE = 'Coef_estim_from_MHE_' + Case + '_Nt_' + Nt + '_Noise_' + Noise +'.txt'
data_MHE2RBD = np.loadtxt(file_name_data_MHE, delimiter=',', skiprows=3)
data_MHE2RBD = data_MHE2RBD[0:len(data_MHE2RBD[:,0])-1,:]

legendas = ['M','Ca','Cad2','Cna','Cypa','Clp','Cma','Cmq']#,'Cnpa']

for i in legendas:
    j = legendas.index(i)
    plt.figure()
    plt.plot(data_eg[:,0],data_eg[:,j],'o-r',label='Egip')
    plt.plot(data_ba[:,0],data_ba[:,j],'*-b',label='Baranw')
    plt.plot(data_NWU[:,0],data_NWU[:,j],'*-m',label='NWU')
    plt.plot(data_MHE2RBD[:,0],data_MHE2RBD[:,j],'*-g',label='MHE2RBD')
    plt.xlabel('Mach')
    plt.ylabel(legendas[j])
    plt.title(legendas[j])
    plt.legend()

plt.show()
