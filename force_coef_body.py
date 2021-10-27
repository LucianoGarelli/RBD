import math as math
import numpy as np

def force_coef_body(mach):

#Spark data
    #Cx0 = 0.221
    #Cx2 = 5.0
    #Cna = 5.83
 #Estimation CFD/RBD
    Cx0 = 0.2387
    Cx2 = 5.933
    Cna = 5.644


    return [Cx0, Cx2, Cna]