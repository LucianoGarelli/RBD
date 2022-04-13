import numpy as np
import matplotlib.pyplot as plt
from csaps import csaps
from scipy.signal import savgol_filter

def data_smoothing(x, y, N, alpha, plot,plot_title):
    '''
    x: interpolated time, unified independet variable
    y: value to be smoothed
    N: cant step to discretiza the independent value
    smooth: eg: smooth=0.999, smoother factor parameter
    plot: "True" or "False", boolean value that shows plot of filtered data
    yhat: return data, smoothed/filtered value
    '''
    data = np.asarray([x,y])
    uniqueValues, indicesList = np.unique(data[:,0], return_index=True)
    data_raw = data[indicesList,:]

    xs = np.linspace(data_raw[0,0], data_raw[-1,0], N)
    #yhat = csaps(xs, y, xs, smooth=0.999)
    #yhat = csaps(x, y, x, smooth=0.999)
    yhat = csaps(x, y, x, smooth=alpha)
    #yhat = savgol_filter(xs, 51, 5)
    if plot:
        fig_size = (12,4)
        fig, axs = plt.subplots(1,1, figsize=fig_size)
        fig.canvas.set_window_title('Datos')
        fig.suptitle('Data Filtered - '+plot_title)
        plt.plot(x, y, label='y')
        plt.plot(x, yhat, color='red',label='yhat')
        plt.legend()
        plt.grid()
        
        #plt.show()

    return yhat
