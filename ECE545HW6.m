clear all
clc
d=0.5;
Ro=48^2/500;
C=0.001;
L=0.0008;
s=tf('s');
R = 0.05;

H = Ro*(1-d) / (s^2*Ro*L*C + s*L + Ro*(1-d)^2);

x = 400 * pi*sqrt(-1);
rip1 = 3*abs(Ro*(1-d) / (x^2*Ro*L*C + x*L + Ro*(1-d)^2))

figure
bode(H)

H = (1-d)*Ro / (s^2*L*C*Ro + s*(C*R*Ro+L) + R + (1-d)^2*Ro);

rip2 = 3*abs((1-d)*Ro / (x^2*L*C*Ro + x*(C*R*Ro+L) + R + (1-d)^2*Ro))

figure
bode(H)

d=0.6424;
Ro=30^2/500;

H= d*Ro/(s^2*Ro*L*C + s*(Ro*R*C+L)+R+Ro);

figure
bode(H)

H= d/(s^2*L*C+s*L/Ro+1);

figure
bode(H,{10^2,10^5})