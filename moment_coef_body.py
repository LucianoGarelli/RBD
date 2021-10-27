import math as math
import numpy as np
from scipy import interpolate

def moment_coef_body(mach):

#Spark data
    #Cma = -12.6
    #Cmq = -196
    #Clp = -2.71
#Estimation CFD/RBD
    Cma = -13.8278
    Cmq = -134.4
    Clp = -3.379
    return [Cma, Cmq, Clp]