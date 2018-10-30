%Leighton Clarke
%EE 537 HW #6, Problem 4
close all
Pgen = [0 0.2 0.4 0.6 0.8 1 .9 .8 .6 .4 .2 0];
Qgen = [-.42 -.44 -.42 -.36 -.22 0 .42 .6 .72 .8 .84 .88];
PUEL = [0 0.2 0.4 0.6 0.8 1];
QUEL = [ -.42 -.4 -.37 -.28 -.18 0];

%Part i: SSSL in P-Q plane
Xd = 1.867;
PF = 0.8;
n = 3600;
Xdp = 0.2;
Xt = 0.32;
Xe = 0.1;
CTR = 1600/5;
VTR = 8000/67;
ZTR = VTR/CTR;
V=1;
j=sqrt(-1);
kk = 180/pi;
S = 29.6e6;
VL = 13.8e3;
Zb = VL^2/S;
Xs = Xt + Xe;
rsssl = V^2/2 * (1/Xs+1/Xd);
Qcsssl = V^2/2 * (1/Xs-1/Xd);
theta = linspace(0,2*pi);

plot(Pgen,Qgen)
hold on
grid on
plot(PUEL,QUEL)
plot(rsssl*cos(theta),rsssl*sin(theta)+Qcsssl)
xlabel('P(pu)')
ylabel('Q(pu)')
title('Generator Curves and SSSL in P-Q Plane')

%Part ii: Convert curves to impedance plane
rsssl2 = 1/2*(Xs+Xd)*Zb/ZTR;
Xsssl = -1/2*(Xd-Xs)*Zb/ZTR;
Zgen = zeros(1,length(Pgen));
for k = 1:length(Pgen)
    Zgen(k) = V^2/(conj(Pgen(k)+j*Qgen(k)))*Zb/ZTR;
end
ZUEL = zeros(1,length(PUEL));
for k = 1:length(PUEL)
    ZUEL(k) = V^2/(conj(PUEL(k)+j*QUEL(k)))*Zb/ZTR;
end
figure
plot(real(Zgen),imag(Zgen))
hold on
grid on
plot(real(ZUEL),imag(ZUEL))
plot(rsssl2*cos(theta),rsssl2*sin(theta)+Xsssl)
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Generator Curves and SSSL in Impedance Plane')

%Part iii: plot inner and outer MHO and directional elements
%Inner MHO:
Dm1 = (1.1*Xd - Xdp/2)*Zb/ZTR;
rm1 = Dm1/2;
Xm1 = -Xdp/2*Zb/ZTR;
Xcm1 = Xm1 - rm1;
%Outer MHO:
Dm2 = (1.1*Xd+Xt+Xe)*Zb/ZTR;
rm2 = Dm2/2;
Xm2 = (Xe+Xt)*Zb/ZTR;
Xcm2 = Xm2 - rm2;
%Directional element:
xdir = linspace(-30,30);
ydir = xdir*tan(-30/kk);

%Part iv: plot all curves together in impedance plane
figure
plot(rm1*cos(theta),rm1*sin(theta)+Xcm1)
hold on
grid on
plot(xdir,ydir)
plot(rm2*cos(theta),rm2*sin(theta)+Xcm2)
plot(real(ZUEL),imag(ZUEL))
plot(rsssl2*cos(theta),rsssl2*sin(theta)+Xsssl)
plot(real(Zgen),imag(Zgen))
ylim([-50,25])
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('All Generator and Relay Characteristics in Impedance Plane')