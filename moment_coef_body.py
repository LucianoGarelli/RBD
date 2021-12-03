import math as math
import numpy as np
from scipy import interpolate

def moment_coef_body(mach):

############################
    if 1:
# Costello-Sahu AIAA 2007-6582
#Spark data
        #Cma = -12.6
        #Cmq = -196
        #Clp = -2.71
#Estimation CFD/RBD
        Cma = -13.8278
        Cmq = -134.4
        Clp = -3.379
##########################
#Stahl-Costello-Sahu AIAA 2009-5715
#Synthetic
    if 0:
        Cma = -8.9
        Cmq = -198
        Clp = -2.6

    return [Cma, Cmq, Clp]