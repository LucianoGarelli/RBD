import numpy as np
import matplotlib.pyplot as plt
from csaps import csaps
from scipy.signal import savgol_filter

#data = np.loadtxt('Resu_ref/Wernert_AIAA2010_7460/beta_BFP.csv', delimiter=',', skiprows=0)
data = np.loadtxt('./Resu_ref/Wernert_AIAA2010_7460/Caso_F01_unificated/Forces_proc.txt', delimiter=',', skiprows=1)
#data = np.loadtxt('./Resu_ref/Wernert_AIAA2010_7460/Caso_F01_unificated/Moments_proc.txt', delimiter=',', skiprows=1)

uniqueValues, indicesList = np.unique(data[:,0], return_index=True)
data_raw = data[indicesList,:]
i=15

xs = np.linspace(data_raw[0,0], data_raw[-1,0], 5000)
ys = csaps(data_raw[:,0], data_raw[:,i], data_raw[:,0], smooth=0.999)
smoothing_result = csaps(data_raw[:,0], data_raw[:,i], xs)

yhat = savgol_filter(data_raw[:,i], 51, 5)


fig_size = (12,4)
fig, axs = plt.subplots(1,1, figsize=fig_size)
fig.canvas.set_window_title('Datos')
fig.suptitle('Datos filtrados')
plt.plot(data_raw[:,0], data_raw[:,i])
plt.plot(data_raw[:,0],yhat, color='green')
plt.plot(data_raw[:,0], ys, color='red')

plt.grid()

plt.show()