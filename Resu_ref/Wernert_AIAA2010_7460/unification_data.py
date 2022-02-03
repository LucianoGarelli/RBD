# script for unification data from F01 - W7640 to the same time scale. BFP model
import numpy as np
import matplotlib.pyplot as plt
import os

# define unification time
t0 = 0
tf = 97
n = 10000
time = np.linspace(t0,tf,n)

# input data
# alpha
data = np.loadtxt('alpha_BFP.csv', delimiter=',', skiprows=1)
t_alpha = data[:,0]
alpha = data[:,1]
# interpolated data
alpha_inter = np.interp(time, t_alpha, alpha)
data = np.zeros(1) # FIXME, como hago para borrar, la variable, xq si leo un data de 1000 y desp uno de 100, los ultimos 900 valores me quedan del data anterior

# input data
# beta
data = np.loadtxt('beta_BFP.csv', delimiter=',', skiprows=1)
t_beta = data[:,0]
beta = data[:,1]
# interpolated data
beta_inter = np.interp(time, t_beta, beta)
data = np.zeros(1)

# input data
# alpha total
data = np.loadtxt('alpha_tot_BFP.csv', delimiter=',', skiprows=1)
t_alpha_tot = data[:,0]
alpha_tot = data[:,1]
# interpolated data
alpha_tot_inter = np.interp(time, t_alpha_tot, alpha_tot)
data = np.zeros(1)

# input data
# Force X tot
data = np.loadtxt('Forces_X_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fxtot = data[:,0]
Fxtot = data[:,1]
# interpolated data
Fxtot_inter = np.interp(time, t_Fxtot, Fxtot)
data = np.zeros(1)

# input data
# Force Gx
data = np.loadtxt('Forces_X_G_BFP.csv', delimiter=',', skiprows=1)
t_Fxg = data[:,0]
Fxg = data[:,1]
# interpolated data
Fxg_inter = np.interp(time, t_Fxg, Fxg)
data = np.zeros(1)

Fx = Fxtot_inter-Fxg_inter # hay que sustraerle el aporte de la gravedad, ya que el modelo de estimacion no lo contempla

# input data
# Force Y tot
data = np.loadtxt('Forces_Y_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fytot = data[:,0]
Fytot = data[:,1]
# interpolated data
Fytot_inter = np.interp(time, t_Fytot, Fytot)
data = np.zeros(1)

# input data
# Force Gy - es zero
#data = np.loadtxt('Forces_Y_G_BFP.csv', delimiter=',', skiprows=1)
#t_Fyg = data[:,0]
#Fyg = data[:,1]
# interpolated data
#Fyg_inter = np.interp(time, t_Fyg, Fyg)
pp =  np.shape(Fytot_inter)
Fyg_inter = np.zeros(pp)
data = np.zeros(1)

Fy = Fytot_inter-Fyg_inter

# input data
# Force Z tot
data = np.loadtxt('Forces_Z_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fztot = data[:,0]
Fztot = data[:,1]
# interpolated data
Fztot_inter = np.interp(time, t_Fztot, Fztot)
data = np.zeros(1)

# input data
# Force Gz
data = np.loadtxt('Forces_Z_G_BFP.csv', delimiter=',', skiprows=1)
t_Fzg = data[:,0]
Fzg = data[:,1]
# interpolated data
Fzg_inter = np.interp(time, t_Fzg, Fzg)
data = np.zeros(1)

Fz = Fztot_inter-Fzg_inter

# Vel NED
data = np.loadtxt('mag_V_BFP.csv', delimiter=',', skiprows=1)
t_vel = data[:,0]
vel = data[:,1]
# interpolated data
vel_inter = np.interp(time, t_vel, vel)
data = np.zeros(1)

# p -roll
# spin rate BFP = spin rate BR
data = np.loadtxt('spin_rate_BR.csv', delimiter=',', skiprows=1)
t_p = data[:,0]
p = data[:,1]
# interpolated data
p_inter = np.interp(time, t_p, p)
data = np.zeros(1)

# q - pith rate
data = np.loadtxt('pitch_rate_q_BFP.csv', delimiter=',', skiprows=1)
t_q = data[:,0]
q = data[:,1]
# interpolated data
q_inter = np.interp(time, t_q, q)
data = np.zeros(1)

# r -yaw rate
data = np.loadtxt('yaw_rate_BFP.csv', delimiter=',', skiprows=1)
t_r = data[:,0]
r = data[:,1]
# interpolated data
r_inter = np.interp(time, t_r, r)
data = np.zeros(1)

# z - altura
data = np.loadtxt('ZE_QFW.csv', delimiter=',', skiprows=1)
t_ze = data[:,0]
ze = data[:,1]
# interpolated data
ze_inter = np.interp(time, t_ze, ze)
data = np.zeros(1)

# Moment X
data = np.loadtxt('Moment_X_tot_BFP.csv', delimiter=',', skiprows=1)
t_mx = data[:,0]
mx = data[:,1]
# interpolated data
mx_inter = np.interp(time, t_mx, mx)
data = np.zeros(1)

# Moment Y
data = np.loadtxt('Moment_Y_tot_BFP.csv', delimiter=',', skiprows=1)
t_my = data[:,0]
my = data[:,1]
# interpolated data
my_inter = np.interp(time, t_my, my)
data = np.zeros(1)

# Moment Z
data = np.loadtxt('Moment_Z_tot_BFP.csv', delimiter=',', skiprows=1)
t_mz = data[:,0]
mz = data[:,1]
# interpolated data
mz_inter = np.interp(time, t_mz, mz)
data = np.zeros(1)

# File to write force - Body Frame
ff = open("./Forces_proc_F01_unificated.txt", "w")
ff.write("# Time,     alpha,     beta,     V_inf,    V_inf, V_inf, V_inf,  p, "
         "      q,         r,         gx,         gy,       gz,      FX,         FY,         FZ       ,       altura ZE\n")
ff.close()

ff = open('./Forces_proc_F01_unificated.txt', 'ab')
#f_force = np.asarray([dt*(k+1), alfa, beta, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], C_body[0],C_body[1],C_body[2]])
f_force = np.asarray([time, alpha_inter, beta_inter, vel_inter, vel_inter, vel_inter, vel_inter, p_inter, q_inter, r_inter, Fxg_inter, Fyg_inter, Fzg_inter, Fx, Fy, Fz, ze_inter])
f_force = f_force.T
#np.savetxt(ff, [f_force], delimiter=", ", fmt='%1.3e')
np.savetxt(ff, f_force, delimiter=", ", fmt='%1.3e')
ff.close()

# File to write moment - Body Frame
fm = open("./Moments_proc_F01_unificated.txt", "w")
fm.write("# Time,     alpha,      beta,     p,      q,      r,    MX,     MY,     MZ \n")
fm.close()

fm = open('./Moments_proc_F01_unificated.txt', 'ab')
# modificar este para momentos wind
m_moment = np.asarray([time, alpha_inter, beta_inter, p_inter,  q_inter, r_inter, mx_inter, my_inter, mz_inter])
m_moment = m_moment.T
np.savetxt(fm, m_moment, delimiter=", ", fmt='%1.3e')
fm.close()



'''
#print(np.shape(data))
plt.plot(t_alpha,alpha,label='alpha BFP')
plt.plot(time,alpha_inter,label='alpha inter')
plt.legend()
#plt.set_xlim([min(time), max(time)])
plt.title('Alpha BFP')
plt.xlabel('Time [s]')
'''
#print(np.shape(data))
plt.plot(t_vel,vel,label='ZE BFP')
plt.plot(time,vel_inter,label='ZE inter')
plt.legend()
#plt.set_xlim([min(time), max(time)])
plt.title('ZE BFP')
plt.xlabel('Time [s]')


print('*-----------------------------------------------*')
print('\nFIN! - OK.\n')

plt.show()
