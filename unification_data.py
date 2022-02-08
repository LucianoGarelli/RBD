# script for unification data from F01 - W7640 to the same time scale. BFP model
import numpy as np
import matplotlib.pyplot as plt
import os
import save_data as sv

# define unification time
t0 = 10
tf = 97
N = 1000
time = np.linspace(t0,tf,N)

# input data
# alpha
data = []
dir_save = './Resu_ref/Wernert_AIAA2010_7460/Caso_F01_unificated_V2/'
dir_load = './Resu_ref/Wernert_AIAA2010_7460/'

data = np.loadtxt(dir_load + 'alpha_BFP.csv', delimiter=',', skiprows=1)
t_alpha = data[:,0]
alpha = -1*data[:,1]
# interpolated data
alpha_inter = np.interp(time, t_alpha, alpha)
alpha_inter = np.deg2rad(alpha_inter)

# input data
# beta
data = []
data = np.loadtxt(dir_load + 'beta_BFP.csv', delimiter=',', skiprows=1)
t_beta = data[:,0]
beta = -1*data[:,1]
# interpolated data
beta_inter = np.interp(time, t_beta, beta)
beta_inter = np.deg2rad(beta_inter)

# input data
# alpha total
data = []
data = np.loadtxt(dir_load + 'alpha_tot_BFP.csv', delimiter=',', skiprows=1)
t_alpha_tot = data[:,0]
alpha_tot = data[:,1]
# interpolated data
alpha_tot_inter = np.interp(time, t_alpha_tot, alpha_tot)

# input data
# Force X tot
data = []
data = np.loadtxt(dir_load +'Forces_X_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fxtot = data[:,0]
Fxtot = data[:,1]
# interpolated data
Fxtot_inter = np.interp(time, t_Fxtot, Fxtot)

# input data
# Force Gx
data = []
data = np.loadtxt(dir_load +'Forces_X_G_BFP.csv', delimiter=',', skiprows=1)
t_Fxg = data[:,0]
Fxg = data[:,1]
# interpolated data
Fxg_inter = np.interp(time, t_Fxg, Fxg)

Fx = Fxtot_inter-Fxg_inter # hay que sustraerle el aporte de la gravedad, ya que el modelo de estimaci'on no lo contempla

# input data
# Force Y tot
data = []
data = np.loadtxt(dir_load +'Forces_Y_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fytot = data[:,0]
Fytot = data[:,1]
# interpolated data
Fytot_inter = np.interp(time, t_Fytot, Fytot)

# input data
# Force Gy - es zero
#data = np.loadtxt(dir_load +'Forces_Y_G_BFP.csv', delimiter=',', skiprows=1)
#t_Fyg = data[:,0]
#Fyg = data[:,1]
# interpolated data
#Fyg_inter = np.interp(time, t_Fyg, Fyg)
pp =  np.shape(Fytot_inter)
Fyg_inter = np.zeros(pp)

Fy = Fytot_inter-Fyg_inter

# input data
# Force Z tot
data = []
data = np.loadtxt(dir_load +'Forces_Z_tot_BFP.csv', delimiter=',', skiprows=1)
t_Fztot = data[:,0]
Fztot = data[:,1]
# interpolated data
Fztot_inter = np.interp(time, t_Fztot, Fztot)

# input data
# Force Gz
data = []
data = np.loadtxt(dir_load +'Forces_Z_G_BFP.csv', delimiter=',', skiprows=1)
t_Fzg = data[:,0]
Fzg = data[:,1]
# interpolated data
Fzg_inter = np.interp(time, t_Fzg, Fzg)

Fz = Fztot_inter-Fzg_inter

# p -roll
# spin rate BFP = spin rate BR
data = []
data = np.loadtxt(dir_load +'spin_rate_BR.csv', delimiter=',', skiprows=1)
t_p = data[:,0]
p = data[:,1]
# interpolated data
p_inter = 2.*np.pi*np.interp(time, t_p, p)
#p_inter = np.pi*np.interp(time, t_p, p)

# q - pith rate
data = []
data = np.loadtxt(dir_load +'pitch_rate_q_BFP.csv', delimiter=',', skiprows=1)
t_q = data[:,0]
q = data[:,1]
# interpolated data
q_inter = 2.*np.pi*np.interp(time, t_q, q)
#q_inter = np.pi*np.interp(time, t_q, q)

# r -yaw rate
data = []
data = np.loadtxt(dir_load +'yaw_rate_BFP.csv', delimiter=',', skiprows=1)
t_r = data[:,0]
r = data[:,1]
# interpolated data
r_inter = 2.*np.pi*np.interp(time, t_r, r)
#r_inter = np.pi*np.interp(time, t_r, r)

# vel "x" en NED
data = []
data = np.loadtxt(dir_load +'UE_BFP.csv', delimiter=',', skiprows=1)
t_vxned = data[:,0]
vxned = data[:,1]
# interpolated data
vxned_inter = np.interp(time, t_vxned, vxned)

# vel "y" en NED
data = []
data = np.loadtxt(dir_load +'VE_BFP.csv', delimiter=',', skiprows=1)
t_vyned = data[:,0]
vyned = data[:,1]
# interpolated data
vyned_inter = np.interp(time, t_vyned, vyned)

# vel "z" en NED
data = []
data = np.loadtxt(dir_load +'we_BFP.csv', delimiter=',', skiprows=1)
t_vzned = data[:,0]
vzned = data[:,1]
# interpolated data
vzned_inter = np.interp(time, t_vzned, vzned)

# Vel NED
data = []
data = np.loadtxt(dir_load +'mag_V_BFP.csv', delimiter=',', skiprows=1)
t_vel = data[:,0]
vel = data[:,1]
# interpolated data
vel_inter = np.interp(time, t_vel, vel)

vt_calc = np.sqrt(3)
vt_calc = np.sqrt(vxned_inter**2+vyned_inter**2+vzned_inter**2)
'''
# comparison mag vel NED, digitalized values vs calculated by components
plt.plot(t_vel,vel,label = 'mag V graf')
plt.plot(time,vt_calc,label = 'mag V calc')
plt.legend()
plt.title('Comparison Mag Vel NED')
plt.xlabel('Time [s]')
plt.ylabel('Vt [m/s]')
'''
# FIXME ! ponemos valores pero en realidad abria que calcularlos
ue_inter = vel_inter
ve_inter = vel_inter
we_inter = vel_inter

'''
# la vel est'a en NED, pero Data.hdf5 guarda en body y desp post-procesa
# vel body
vel_body = np.empty((N+1,3))
for k, xk in enumerate(x):
     quat = utiles.Quaternion(xk[6], xk[7:10])
     quat_trans = quat.T
     vel_body[k] = quat_trans.rotate_vector(vel_body)
'''


# x - recorrido longitud
data = []
data = np.loadtxt(dir_load +'XE_QFW.csv', delimiter=',', skiprows=1)
t_xe = data[:,0]
xe = data[:,1]
# interpolated data
xe_inter = 1000.*np.interp(time, t_xe, xe)

# y - recorrido  lateral
data = []
data = np.loadtxt(dir_load +'YE_QFW.csv', delimiter=',', skiprows=1)
t_ye = data[:,0]
ye = data[:,1]
# interpolated data
#ye_inter = -1.*np.interp(time, t_ye, ye)
ye_inter = np.interp(time, t_ye, ye)

# z - altura
data = []
data = np.loadtxt(dir_load +'ZE_QFW.csv', delimiter=',', skiprows=1)
t_ze = data[:,0]
ze = data[:,1]
# interpolated data
ze_inter = 1000*np.interp(time, t_ze, ze)

# Moment X
data = []
data = np.loadtxt(dir_load +'Moment_X_tot_BFP.csv', delimiter=',', skiprows=1)
t_mx = data[:,0]
mx = data[:,1]
# interpolated data
mx_inter = np.interp(time, t_mx, mx)

# Moment Y
data = []
data = np.loadtxt(dir_load +'Moment_Y_tot_BFP.csv', delimiter=',', skiprows=1)
t_my = data[:,0]
my = data[:,1]
# interpolated data
my_inter = np.interp(time, t_my, my)

# Moment Z
data = []
data = np.loadtxt(dir_load + 'Moment_Z_tot_BFP.csv', delimiter=',', skiprows=1)
t_mz = data[:,0]
mz = data[:,1]
# interpolated data
mz_inter = np.interp(time, t_mz, mz)

CHECK_FOLDER = os.path.isdir(dir_save)
# If folder doesn't exist, then create it.
if not CHECK_FOLDER:
    os.makedirs(dir_save)
    print("Directorio creado: ", dir_save)
else:
    print(dir_save, "Directorio existente.")

# File to write force - Body Frame

ff = open(dir_save + 'Forces_proc.txt', 'w')
ff.write("# Time,     alpha,     beta,     V_inf,    V_inf, V_inf, V_inf,  p, "
         "      q,         r,         gx,         gy,       gz,      FX,         FY,         FZ       \n")
ff.close()

ff = open( dir_save +'Forces_proc.txt', 'ab')
#f_force = np.asarray([dt*(k+1), alfa, beta, vt, x[3], x[4], x[5], x[10], x[11], x[12],g_body[0],g_body[1],g_body[2], C_body[0],C_body[1],C_body[2]])
f_force = np.asarray([time, alpha_inter, beta_inter, vel_inter, ue_inter, ve_inter, we_inter, p_inter, q_inter, r_inter, Fxg_inter, Fyg_inter, Fzg_inter, Fx, Fy, Fz])
f_force = f_force.T
np.savetxt(ff, f_force, delimiter=", ", fmt='%1.3e')
ff.close()

# File to write moment - Body Frame
fm = open(dir_save + 'Moments_proc.txt', 'w')
fm.write("# Time,     alpha,      beta,     p,      q,      r,    MX,     MY,     MZ \n")
fm.close()

fm = open(dir_save + 'Moments_proc.txt', 'ab')
# modificar este para momentos wind
m_moment = np.asarray([time, alpha_inter, beta_inter, p_inter,  q_inter, r_inter, mx_inter, my_inter, mz_inter])
m_moment = m_moment.T
np.savetxt(fm, m_moment, delimiter=", ", fmt='%1.3e')
fm.close()

N = np.size(p_inter)

# creamos le vector de estados, para usar save_data
x = np.zeros((N,15))
#x[:,0:3] = [ze_inter ze_inter ze_inter]
x[:,0] = xe_inter
x[:,1] = ye_inter
x[:,2] = ze_inter
# aca empieza a rellenar CHECK
x[:,3] = vt_calc #vel_body[0]
#x[:,4] = vt_calc#vel_body[1]
#x[:,5] = vt_calc#vel_body[2]
'''
x[:,6] = ze_inter
x[:,7] = ze_inter
x[:,8] = ze_inter
x[:,9] = ze_inter
'''
#x[3] = V*math.cos(alfa)*math.cos(beta) # velocidad inicial en body
#x[4] = V*math.sin(beta)
#x[5] = V*math.sin(alfa)*math.cos(beta)
x[:,10] = p_inter # rotacion inicial en body
x[:,11] = q_inter
x[:,12] = r_inter
x[:,13] = alpha_inter
x[:,14] = beta_inter

print(x[1,:])
'''
#print(np.shape(data))
plt.plot(t_alpha,alpha,label='alpha BFP')
plt.plot(time,alpha_inter,label='alpha inter')
plt.legend()
#plt.set_xlim([min(time), max(time)])
plt.title('Alpha BFP')
plt.xlabel('Time [s]')

#print(np.shape(data))
plt.plot(t_Fztot,Fztot,label='Fz tot BFP')
plt.plot(time,Fztot_inter,label='Fz tot inter')
plt.plot(time,Fz,label='Fz s/G inter')
plt.plot(time,Fzg_inter,label='Fz G inter')
plt.legend()
#plt.set_xlim([min(time), max(time)])
plt.title('FZ BFP')
plt.xlabel('Time [s]')
'''

plt.plot(t_mz,mz,label='Mz tot BFP')
#plt.plot(time,Fztot_inter,label='Fz tot inter')
#plt.plot(time,Fz,label='Fz s/G inter')
#plt.plot(time,Fzg_inter,label='Fz G inter')
plt.legend()
#plt.set_xlim([min(time), max(time)])
plt.title('MZ BFP')
plt.xlabel('Time [s]')

#sv.save_data(N, time, x, Ixx, Iyy, Izz)
sv.save_data(N-1, time, x, 0, 0, 0,dir_save)

print('*-----------------------------------------------*')
print('\nFIN! - OK.\n')

plt.show()
