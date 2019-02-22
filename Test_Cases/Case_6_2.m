# Example 5.2 Modern Exterior Ballistics

Cd = 0.205
theta = 0.2099
V = 860
t = 0.756
rho = 0.0765
S = 0.001114
m = 0.03286
g = 32.1740
Vx0 = V*cosd(theta)
Y0 = -0.05
  
k1 = (rho*S*Cd)/(2*m)

x = log(t*Vx0*k1+1)/k1 #Carlucci
Vxy = Vx0*exp(-k1*x) 
x1 = t*log(Vx0/Vxy)*Vx0/(Vx0/Vxy - 1) #Macoy 5.43


Vz = Vxy*[tand(theta) - (g*t)*(1 + Vx0*k1*t/2)/(Vx0)]
a = (t^2*g/2)*(0.5 + (Vx0/Vxy - 1)^(-1) - (Vx0/Vxy - 1)^(-2)*log(Vx0/Vxy) );

Y = Y0 + x*tand(theta) - a
