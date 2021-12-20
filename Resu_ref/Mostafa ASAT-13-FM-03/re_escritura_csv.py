import numpy as np
import matplotlib.pyplot as plt

with open ("./Vel_mag.csv") as f_input:
  text = [l.replace(",", ".") for l in f_input]

data = np.loadtxt(text, delimiter=';')

np.savetxt('Vel_mag.csv',data,delimiter=";")
