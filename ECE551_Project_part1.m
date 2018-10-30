clear all
close all

%define parameters given in circuit description

Vg = 12;
Vo = 24;
r = 28;
c1 = 2*10^-6;
c2 = 2*10^-5;
d = Vo/(Vg+Vo);
l1 = 5*10^-4;
r1 = 0.01;
l2 = 7.5*10^-3;
r2 = 0.01;
m = -1.5*10^-3;
s2 = l1*l2-m^2;

%Define ss variables for switch on(1) and off(2)

A1=[-1/(r*c2), 0, 1/c2, 0;
    0, 0, -1/c1, 0;
    -l1/s2, l1/s2, -l1*r2/s2, m*r1/s2;
    m/s2, -m/s2, m*r2/s2, -l2*r1/s2];

B1=[0;0;-m/s2;l2/s2];
C1=[1,0,0,0];
D1=0;

A2=[-1/(r*c2), 0, 1/c2, 0;
    0, 0, 0, 1/c1; %This 1/c1 term is negative in book!!!
    -l1/s2, m/s2, -l1*r2/s2, m*r1/s2;
    m/s2, -l2/s2, m*r2/s2, -l2*r1/s2];
B2=B1;
C2=C1;
D2=D1;

%Use duty cycle to create ss variables for averaged model
d1 = 1-d;
A = d*A1+d1*A2;
B = d*B1+d1*B2;
C = d*C1+d1*C2;

%Create Bd variable for control input
X = -A^-1*B*Vg;
Bd = (A1-A2)*X;