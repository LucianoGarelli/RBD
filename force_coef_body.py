import math as math
import numpy as np

def force_coef_body(mach):

############################
    if 1:
# Costello-Sahu AIAA 2007-6582
#Spark data
        #Cx0 = 0.221
        #Cx2 = 5.0
        #Cna = 5.83
 #Estimation CFD/RBD
        Cx0 = 0.2387
        Cx2 = 5.933
        Cna = 5.644
##########################
#Stahl-Costello-Sahu AIAA 2009-5715
#Synthetic
    if 0:
        Cx0 = 0.23
        Cx2 = 5
        Cna = 5.9


    return [Cx0, Cx2, Cna]