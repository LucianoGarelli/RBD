# Ecuacione Carlucci Pag. 207
# Flat-fire theta<5.7 deg

Cd = 0.3
theta = 5
V = 1000
t = 6
rho = 1.225
S = pi*0.02^2
m = 0.1
g = 9.81
Vx0 = V*cosd(theta)
Y0 = 10
  
k1 = (rho*S*Cd)/(2*m)
x = log(t*Vx0*k1+1)/k1 #Carlucci

Vxy = Vx0*exp(-k1*x)
Vxy2 = Vx0/(1+Vx0*k1*t)

x1 = t*log(Vx0/Vxy)*Vx0/(Vx0/Vxy - 1) #Mccoy 5.43
Vz = Vxy*[tand(theta) - (g*t)*(1 + Vx0*k1*t/2)/(Vx0)]
a = (t^2*g/2)*(0.5 + (Vx0/Vxy - 1)^(-1) - (Vx0/Vxy - 1)^(-2)*log(Vx0/Vxy) );

Y = Y0 + x*tand(theta) - a



# octave(x=1167.0,Vx=67.303, Vz=25.53, Y=13.515)
# python(x=1164.6,Vx=66.532, Vz=25.356, Y=13.475)
