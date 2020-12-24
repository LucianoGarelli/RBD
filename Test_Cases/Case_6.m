# Ecuacione Carlucci Pag. 207
# Flat-fire theta<5.7 deg
clear all;close all

Cd = 0.1
theta = 45
V = 100
t = 12.643
rho = 1.225
diam = 0.04
S = pi/4*diam^2
m = 0.1
g = 9.81
Vx0 = V*cosd(theta)
Y0 = 0
  
k1 = (rho*S*Cd)/(2*m)
x = log(t*Vx0*k1+1)/k1 #Carlucci

% velocidad en tiempo final - Vel comp "x"
Vxy = Vx0*exp(-k1*x)
Vxy2 = Vx0/(1+Vx0*k1*t) % idem anterior, otra formula

x1 = t*log(Vx0/Vxy)*Vx0/(Vx0/Vxy - 1) #Mccoy 5.43
% velocidad en tiempo final - Vel comp "z"
Vz = Vxy*[tand(theta) - (g*t)*(1 + Vx0*k1*t/2)/(Vx0)]

a = (t^2*g/2)*(0.5 + (Vx0/Vxy - 1)^(-1) - (Vx0/Vxy - 1)^(-2)*log(Vx0/Vxy) );

Y = Y0 + x*tand(theta) - a

if 0 % parte symbolic 
  syms Cd rho S m t Vx0
  k1 = (rho*S*Cd)/(2*m)
  xs = log(t*Vx0*k1+1)/k1 %Carlucci

  a = (t^2*g/2)*(0.5 + (Vx0/Vxy - 1)^(-1) - (Vx0/Vxy - 1)^(-2)*log(Vx0/Vxy) );
  Y = Y0 + x*tand(theta) - a

  xs_subs=subs(xs,{m,Cd,S,Vx0,rho,t},{0.1,0.1,0.0012566,70.711,1.225,12.6428})
end

# octave(x=1167.0,Vx=67.303, Vz=25.53, Y=13.515)
# python(x=1164.6,Vx=66.532, Vz=25.356, Y=13.475)
