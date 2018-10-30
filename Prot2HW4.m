%Leighton Clarke
%Advanced Protections, Homework 4-5
%May 7, 2018

close all
%%Problem 1, part a
kk = 180/pi;
a = 1*exp(j*120/kk);
j = sqrt(-1);
CTR = 800/5;
VTR = 500000/sqrt(3)/ 67;
ZTR = VTR/CTR;

Zs1 = 1 + j*10;
Zs2 = Zs1;
Zs0 = 2 + j*30;
ZR1 = j*20;
ZR2 = ZR1;
ZR0 = j*10;
zL1 = 0.073 + j*0.8;
zL2 = zL1;
zL0 = 0.1 + j*2.6;
ZL1B23 = 100*zL1;
ZL1B12 = 50*zL1;
ZL0B23 = zL0*100;
ZL0B12 = zL0*50;
d1 = 20;
d2 = 10;
m1 = d1/100;
m2 = d2/50;
Zr = ZL1B23/ZTR*0.8;
Vbase = 500000;
Es = 1.002*exp(j*12.5/kk)*Vbase/sqrt(3);
ER = 0.994*Vbase/sqrt(3);
Ipf = (Es-ER)/(Zs1+ZL1B12+ZL1B23+ZR1); %Pre-fault current from B1 to B3

%Case 1: a-g fault at F1 with Rf = 0 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f1agr0 = csvread('f1_ag_r0.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f1agr0,CTR,VTR);
fprintf('\nProblem 1, part a:\n')
fprintf('\n-----Case 1: a-g fault at F1 with Rf = 0 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
Zx1 = m1*ZL1B23 + Zs1 + ZL1B12;
Zy1 = (1-m1)*ZL1B23 + ZR1;
Z1 = Zx1*Zy1/(Zx1+Zy1);
Zx0 = m1*ZL0B23 + Zs0 + ZL0B12;
Zy0 = (1-m1)*ZL0B23 + ZR0;
Z0 = Zx0*Zy0/(Zx0+Zy0);
C1 = Zy1/(Zx1+Zy1);
C0 = Zy0/(Zx0+Zy0);
Eth = C1*Es+(1-C1)*ER;
If1 = Eth/(2*Z1+Z0);
alpha = Ipf/If1;
I1 = (C1+alpha)*If1;
I2 = C1*If1;
I0 = C0*If1;
IA = I1+I2+I0;
IB = a^2*I1+a*I2+I0;
IC = a*I1 + a^2*I2 +I0;
V1 = If1*((C1+alpha)*(m1*ZL1B23)+Z1+Z0);
V2 = If1*(C1*(m1*ZL1B23)-Z1);
V0 = If1*(C0*(m1*ZL0B23)-Z0);
%VA = If1*((2*C1+alpha)*(m1*ZL1B23+ZL1B12)+C0*(m1*ZL0B23+ZL0B12));
VA = V1+V2+V0;
VB = a^2*V1+a*V2+V0;
VC = a*V1+a^2*V2+V0;
Va1 = VA/VTR;
Vb1 = VB/VTR;
Vc1 = VC/VTR;
Ia1 = IA/CTR;
Ib1 = IB/CTR;
Ic1 = IC/CTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va1), angle(Va1)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb1), angle(Vb1)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc1), angle(Vc1)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia1), angle(Ia1)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib1), angle(Ib1)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic1), angle(Ic1)*kk)

%Case 2: a-g fault at F1 with Rf = 2 ohms
%Use digital relay algorithm to obtain voltage and current phasors
f1agr2 = csvread('f1_ag_r2.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f1agr2,CTR,VTR);

fprintf('\n-----Case 2: a-g fault at F1 with Rf = 2 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
Rf = 2;
If1 = Eth/(2*Z1+Z0+3*Rf);
alpha = Ipf/If1;
I1 = (C1+alpha)*If1;
I2 = C1*If1;
I0 = C0*If1;
IA = I1+I2+I0;
IB = a^2*I1+a*I2+I0;
IC = a*I1 + a^2*I2 +I0;
V1 = If1*((C1+alpha)*(m1*ZL1B23)+Z1+Z0+3*Rf);
V2 = If1*(C1*(m1*ZL1B23)-Z1);
V0 = If1*(C0*(m1*ZL0B23)-Z0);
VA = V1+V2+V0;
VB = a^2*V1+a*V2+V0;
VC = a*V1+a^2*V2+V0;
Va2 = VA/VTR;
Vb2 = VB/VTR;
Vc2 = VC/VTR;
Ia2 = IA/CTR;
Ib2 = IB/CTR;
Ic2 = IC/CTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va2), angle(Va2)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb2), angle(Vb2)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc2), angle(Vc2)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia2), angle(Ia2)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib2), angle(Ib2)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic2), angle(Ic2)*kk)

%Case 3: b-c fault at F1 with Rf = 0 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f1bcr0 = csvread('f1_bc_r0.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f1bcr0,CTR,VTR);

fprintf('\n-----Case 3: b-c fault at F1 with Rf = 0 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
If1 = Eth/(2*Z1);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = -C1*If1;
IA = I1+I2;
IB = a^2*I1 + a*I2;
IC = a*I1 + a^2*I2;
V1 = If1*((C1+alpha)*m1*ZL1B23+Z1);
V2 = If1*(-C1*m1*ZL1B23+Z1);
VA = V1+V2;
VB = a^2*V1 + a*V2;
VC = a*V1 + a^2*V2;
Va3 = VA/VTR;
Vb3 = VB/VTR;
Vc3 = VC/VTR;
Ia3 = IA/CTR;
Ib3 = IB/CTR;
Ic3 = IC/CTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va3), angle(Va3)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb3), angle(Vb3)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc3), angle(Vc3)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia3), angle(Ia3)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib3), angle(Ib3)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic3), angle(Ic3)*kk)

%Case 4: b-c fault at F1 with Rf = 2 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f1bcr2 = csvread('f1_bc_r2.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f1bcr2,CTR,VTR);

fprintf('\n-----Case 4: b-c fault at F1 with Rf = 2 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
If1 = Eth/(2*Z1+Rf);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = -C1*If1;
IA = I1+I2;
IB = a^2*I1+a*I2;
IC = a*I1+a^2*I2;
V1 = If1*((C1+alpha)*m1*ZL1B23+Z1+Rf);
V2 = If1*(-C1*m1*ZL1B23+Z1);
VA = V1+V2;
VB = a^2*V1+a*V2;
VC = a*V1+a^2*V2;
Va4 = VA/VTR;
Vb4 = VB/VTR;
Vc4 = VC/VTR;
Ia4 = IA/CTR;
Ib4 = IB/CTR;
Ic4 = IC/CTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va4), angle(Va4)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb4), angle(Vb4)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc4), angle(Vc4)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia4), angle(Ia4)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib4), angle(Ib4)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic4), angle(Ic4)*kk)

%Case 5: a-g fault at F2 with Rf = 0 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f2agr0 = csvread('f2_ag_r0.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f2agr0,CTR,VTR);

fprintf('\n-----Case 5: a-g fault at F2 with Rf = 0 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
%Here I "reverse" the system
ER = Es;
Es = Vbase/sqrt(3)*.994;
Zs1eq = ZR1+ZL1B23; %Equivalent source impedance at relay location (includes line B23)
Zs0eq = ZR0+ZL0B23;
Ipf = (Es-ER)/(Zs1+ZL1B12+ZL1B23+ZR1); %Pre-fault current from B3 to B1
Zx1 = Zs1eq + m2*ZL1B12;
Zy1 = (1-m2)*ZL1B12+Zs1;
Zx0 = Zs0eq + m2*ZL0B12;
Zy0 = (1-m2)*ZL0B12 + Zs0;
Z1 = Zx1*Zy1/(Zx1+Zy1);
C1 = Zy1/(Zx1+Zy1);
Z0 = Zx0*Zy0/(Zx0+Zy0);
C0 = Zy0/(Zx0+Zy0);

Eth = (1-C1)*ER+C1*Es;
If1 = Eth/(2*Z1+Z0);
alpha = Ipf/If1;
I1 = If1*(C1 + alpha);
I2 = If1*C1;
I0 = If1*C0;
Ia5 = (I1+I2+I0)/CTR;
Ib5 = (a^2*I1+a*I2+I0)/CTR;
Ic5 = (a*I1+a^2*I2+I0)/CTR;
V1 = If1*((C1+alpha)*m2*ZL1B12+Z1+Z0);
V2 = If1*(C1*m2*ZL1B12-Z1);
V0 = If1*(C0*m2*ZL0B12-Z0);
Va5 = (V1+V2+V0)/VTR;
Vb5 = (a^2*V1+a*V2+V0)/VTR;
Vc5 = (a*V1+a^2*V2+V0)/VTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va5), angle(Va5)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb5), angle(Vb5)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc5), angle(Vc5)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia5), angle(Ia5)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib5), angle(Ib5)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic5), angle(Ic5)*kk)

%Case 6: a-g fault at F2 with Rf = 2 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f2agr2 = csvread('f2_ag_r2.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f2agr2,CTR,VTR);

fprintf('\n-----Case 6: a-g fault at F2 with Rf = 2 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
If1 = Eth/(2*Z1+Z0+3*Rf);
alpha = Ipf/If1;
I1 = If1*(C1 + alpha);
I2 = If1*C1;
I0 = If1*C0;
Ia6 = (I1+I2+I0)/CTR;
Ib6 = (a^2*I1+a*I2+I0)/CTR;
Ic6 = (a*I1+a^2*I2+I0)/CTR;
V1 = If1*((C1+alpha)*m2*ZL1B12+Z1+Z0+3*Rf);
V2 = If1*(C1*m2*ZL1B12-Z1);
V0 = If1*(C0*m2*ZL0B12-Z0);
Va6 = (V1+V2+V0)/VTR;
Vb6 = (a^2*V1+a*V2+V0)/VTR;
Vc6 = (a*V1+a^2*V2+V0)/VTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va6), angle(Va6)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb6), angle(Vb6)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc6), angle(Vc6)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia6), angle(Ia6)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib6), angle(Ib6)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic6), angle(Ic6)*kk)

%Case 7: b-c fault at F2 with Rf = 0 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f2bcr0 = csvread('f2_bc_r0.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f2bcr0,CTR,VTR);

fprintf('\n-----Case 7: b-c fault at F2 with Rf = 0 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
If1 = Eth/(2*Z1);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = If1*-C1;
Ia7 = (I1+I2)/CTR;
Ib7 = (a^2*I1+a*I2)/CTR;
Ic7 = (a*I1+a^2*I2)/CTR;
V1 = If1*((C1+alpha)*m2*ZL1B12+Z1);
V2 = If1*(-C1*m2*ZL1B12+Z1);
Va7 = (V1+V2)/VTR;
Vb7 = (a^2*V1+a*V2)/VTR;
Vc7 = (a*V1+a^2*V2)/VTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va7), angle(Va7)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb7), angle(Vb7)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc7), angle(Vc7)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia7), angle(Ia7)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib7), angle(Ib7)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic7), angle(Ic7)*kk)

%Case 8: b-c fault at F2 with Rf = 2 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f2bcr2 = csvread('f2_bc_r2.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f2bcr2,CTR,VTR);

fprintf('\n-----Case 8: b-c fault at F2 with Rf = 2 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
IF1 = Eth/(2*Z1+Rf);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = -C1*If1;
Ia8 = (I1+I2)/CTR;
Ib8 = (a^2*I1+a*I2)/CTR;
Ic8 = (a*I1+a^2*I2)/CTR;
V1 = If1*((C1+alpha)*m2*ZL1B12+Z1+Rf);
V2 = If1*(-C1*m2*ZL1B12+Z1);
Va8 = (V1+V2)/VTR;
Vb8 = (a^2*V1+a*V2)/VTR;
Vc8 = (a*V1+a^2*V2)/VTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va8), angle(Va8)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb8), angle(Vb8)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc8), angle(Vc8)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia8), angle(Ia8)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib8), angle(Ib8)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic8), angle(Ic8)*kk)

%Case 9: b-c-g fault at F1 with Rf = 0 ohms

%Use digital relay algorithm to obtain voltage and current phasors
f2bcr2 = csvread('f1_bcg_r0.csv');
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f2bcr2,CTR,VTR);

fprintf('\n-----Case 9: b-c-g fault at F1 with Rf = 0 ohms-----\n')
fprintf('Values from digital relay algorithm:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Icp(end)), angle(Icp(end))*kk)

%Calculate phasors using sequence components
Es = 1.002*exp(j*12.5/kk)*Vbase/sqrt(3);
ER = 0.994*Vbase/sqrt(3);
Zx1 = m1*ZL1B23 + Zs1 + ZL1B12;
Zy1 = (1-m1)*ZL1B23 + ZR1;
Z1 = Zx1*Zy1/(Zx1+Zy1);
Zx0 = m1*ZL0B23 + Zs0 + ZL0B12;
Zy0 = (1-m1)*ZL0B23 + ZR0;
Z0 = Zx0*Zy0/(Zx0+Zy0);
C1 = Zy1/(Zx1+Zy1);
C0 = Zy0/(Zx0+Zy0);
Eth = C1*Es+(1-C1)*ER;
Ipf = (Es-ER)/(Zs1+ZL1B12+ZL1B23+ZR1);
If1 = Eth/(Z1+Z1*Z0/(Z1+Z0));
If2 = If1*Z0/(Z1+Z0);
If0 = If1*Z1/(Z1+Z0);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = If2*-C1;
I0 = If0*-C0;
Ia9 = (I1+I2+I0)/CTR;
Ib9 = (a^2*I1+a*I2+I0)/CTR;
Ic9 = (a*I1+a^2*I2+I0)/CTR;
V1 = I1*m1*ZL1B23+If2*Z1;
V2 = -I2*(Zs1+ZL1B12);
V0 = -I0*(Zs0+ZL0B12);
Va9 = (V1+V2+V0)/VTR;
Vb9 = (a^2*V1+a*V2+V0)/VTR;
Vc9 = (a*V1+a^2*V2+V0)/VTR;

fprintf('\nValues from calculations:\n')
fprintf('Va: %f V, angle %f degrees\n', abs(Va9), angle(Va9)*kk)
fprintf('Vb: %f V, angle %f degrees\n', abs(Vb9), angle(Vb9)*kk)
fprintf('Vc: %f V, angle %f degrees\n', abs(Vc9), angle(Vc9)*kk)
fprintf('Ia: %f A, angle %f degrees\n', abs(Ia9), angle(Ia9)*kk)
fprintf('Ib: %f A, angle %f degrees\n', abs(Ib9), angle(Ib9)*kk)
fprintf('Ic: %f A, angle %f degrees\n', abs(Ic9), angle(Ic9)*kk)

%%Problem 1, part b
VB2pf = Es - Ipf*(Zs1+ZL1B12); %Pre-fault voltage at bus B2
SB2pf = 3*VB2pf*conj(Ipf);
fprintf('\nProblem 1, part b:\n')
fprintf('Prefault current: %f A, angle %f degrees\n', abs(Ipf), angle(Ipf)*kk)
fprintf('Prefault power (S) at bus 2: %f MVA, angle %f degrees\n', abs(SB2pf)/10^6,angle(SB2pf)*kk)

%%Problem 1, part c
fprintf('\nProblem 1, part c:\n')
fprintf('See attached document for detailed analysis of results.\n')

%%Problem 1, part d
%Calculate fault 1 distance from data
[Vap,Vbp,Vcp,Iap,Ibp,Icp] = DigiRelay(f1agr0,CTR,VTR);
k0 = (zL0-zL1)/(3*zL1);
Ires = Iap+Ibp+Icp;
Zag = Vap./(Iap+k0*Ires)*ZTR;
m = imag(Zag)./imag(ZL1B23);
fprintf('\nProblem 1, part d: F1 distance calculated from relay data\n')
fprintf('F1 distance using case 1 data: %f miles\n', m(end)*100)
fprintf('See attached plot of m vs. time.\n')
%Plot m versus time
t0 = f1agr0(:,1);
t = linspace(0,t0(end),length(m));
plot(t,m)
xlabel('Time (s)')
ylabel('m')
title('Problem 1, part d: m vs time')

%%Problem 2
fprintf('\nProblem 2:\n')
fprintf('See attached diagram for relay labels used in display.\n')
VTR = 115000/sqrt(3)/67;
CTR = 600/5;
ZTR = VTR/CTR;

Z1_12 = 2.352+j*10.716;
Z0_12 = 3+j*30;
Z1_13 = 0.135+j*1.438;
Z0_13 = 0.5+j*5;
Z1_34 = 0.098+j*1.039;
Z0_34 = 0.2+j*3.5;
Z1_45 = 0.481+j*3.747;
Z0_45 = 1+j*10;
Z1_25 = 1.493+j*6.807;
Z0_25 = 2+j*17.5;
Z1_14 = 0.238+j*2.477;
Z0_14 = 0.8+j*8;

k0_12 = (Z0_12-Z1_12)/(Z1_12*3);
k0_13 = (Z0_13-Z1_13)/(Z1_13*3);
k0_34 = (Z0_34-Z1_34)/(Z1_34*3);
k0_45 = (Z0_45-Z1_45)/(Z1_45*3);
k0_25 = (Z0_25-Z1_25)/(Z1_25*3);
k0_14 = (Z0_14-Z1_14)/(Z1_14*3);

ZrAB1 = (0.8*Z1_12)/ZTR;
ZrAB2 = (Z1_12+0.2*Z1_25)/ZTR;
ZrCD1 = ZrAB1;
ZrCD2 = (Z1_12+0.2*Z1_13)/ZTR;
ZrE1 = (0.8*Z1_13)/ZTR;
ZrE2 = (Z1_13+0.2*Z1_34)/ZTR;
ZrF1 = (0.8*Z1_14)/ZTR;
ZrF2 = (Z1_14+0.2*Z1_34)/ZTR;
ZrG1 = ZrE1;
ZrG2 = (Z1_13+0.2*Z1_14)/ZTR;
ZrH1 = (0.8*Z1_34)/ZTR;
ZrH2 = (Z1_34+0.2*Z1_45)/ZTR;
ZrI1 = ZrH1;
ZrI2 = (Z1_34+0.2*Z1_13)/ZTR;
ZrJ1 = ZrF1;
ZrJ2 = (Z1_14+0.2*Z1_13)/ZTR;
ZrK1 = (0.8*Z1_45)/ZTR;
ZrK2 = (Z1_45+0.2*Z1_25)/ZTR;
ZrL1 = ZrK1;
ZrL2 = (Z1_45+0.2*Z1_34)/ZTR;
ZrM1 = (0.8*Z1_25)/ZTR;
ZrM2 = (Z1_25+0.2*Z1_12)/ZTR;
ZrN1 = ZrM1;
ZrN2 = (Z1_25+0.2*Z1_45)/ZTR;

fprintf('All impedances are in secondary ohms\n')
fprintf('Relay A:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrAB1),angle(ZrAB1)*kk,abs(ZrAB2),angle(ZrAB2)*kk)
fprintf('Relay B:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrAB1),angle(ZrAB1)*kk,abs(ZrAB2),angle(ZrAB2)*kk)
fprintf('Relay C:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrCD1),angle(ZrCD1)*kk,abs(ZrCD2),angle(ZrCD2)*kk)
fprintf('Relay D:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrCD1),angle(ZrCD1)*kk,abs(ZrCD2),angle(ZrCD2)*kk)
fprintf('Relay E:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_13,imag(k0_13),abs(ZrE1),angle(ZrE1)*kk,abs(ZrE2),angle(ZrE2)*kk)
fprintf('Relay F:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_14,imag(k0_14),abs(ZrF1),angle(ZrF1)*kk,abs(ZrF2),angle(ZrF2)*kk)
fprintf('Relay G:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_13,imag(k0_13),abs(ZrG1),angle(ZrG1)*kk,abs(ZrG2),angle(ZrG2)*kk)
fprintf('Relay H:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_34,imag(k0_34),abs(ZrH1),angle(ZrH1)*kk,abs(ZrH2),angle(ZrH2)*kk)
fprintf('Relay I:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_34,imag(k0_34),abs(ZrI1),angle(ZrI1)*kk,abs(ZrI2),angle(ZrI2)*kk)
fprintf('Relay J:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_14,imag(k0_14),abs(ZrJ1),angle(ZrJ1)*kk,abs(ZrJ2),angle(ZrJ2)*kk)
fprintf('Relay K:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_45,imag(k0_45),abs(ZrK1),angle(ZrK1)*kk,abs(ZrK2),angle(ZrK2)*kk)
fprintf('Relay L:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_45,imag(k0_45),abs(ZrL1),angle(ZrL1)*kk,abs(ZrL2),angle(ZrL2)*kk)
fprintf('Relay M:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_25,imag(k0_25),abs(ZrM1),angle(ZrM1)*kk,abs(ZrM2),angle(ZrM2)*kk)
fprintf('Relay N:\nk0 = %f+j%f, Zr1 = %f ohms (angle %f degrees) Zr2 = %f ohms (angle %f degrees)\n', k0_25,imag(k0_25),abs(ZrN1),angle(ZrN1)*kk,abs(ZrN2),angle(ZrN2)*kk)
fprintf('\nFor all relays, the time delay for zone 2 is 0.3 seconds, and for zone 1 is instantaneous.\n')

%%Problem 3
%Part a: rewrite filtering algorithm for differential relay data
fprintf('\nProblem 3, part a:\n')
fprintf('See the function "DigiDiffRelay" - rewritten for differential relay data\n')
%Part b: filter data, calculate alpha, and plot characteristic with alpha to determine operation
fprintf('\nProblem 3, part b:\n')
fprintf('See figures: Relays operate if final alpha value lies outside restraining region.\n')
fprintf('Case 1: "A" relay trips\n')
fprintf('Case 2: "A" relay trips\n')
fprintf('Case 3: "B" and "C" relays trip\n')
fprintf('Case 4: "B" and "C" relays trip\n')
fprintf('Case 5: No relays trip\n')
fprintf('Case 6: No relays trip\n')
fprintf('Case 7: No relays trip\n')
fprintf('Case 8: No relays trip\n')
fprintf('Case 9: "B" and "C" relays trip\n')
r1 = 1/6;
r2 = 6;
B = 98;
CTR = 1200/5;

%Case 1
f1agxx = csvread('f1agxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f1agxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #1')

%Case 2
f1r2agxx = csvread('f1r2agxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f1r2agxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #2')

%Case 3
f1bcxx = csvread('f1bcxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f1bcxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #3')

%Case 4
f1r2bcxx = csvread('f1r2bcxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f1r2bcxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #4')

%Case 5
f2agxx = csvread('f2agxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f2agxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #5')

%Case 6
f2r2agxx = csvread('f2r2agxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f2r2agxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #6')

%Case 7
f2bcx = csvread('f2bcx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f2bcx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #7')

%Case 8
f2r2bcxx = csvread('f2r2bcxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f2r2bcxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #8')

%Case 9
f1bcgxx = csvread('f1bcgxx.csv');
[ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(f1bcgxx,CTR);
alpha1a = ISA./IRA;
alpha1b = ISB./IRB;
alpha1c = ISC./IRC;

figure
Rtheta = (-B/kk):0.01:(B/kk);
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha1a),imag(alpha1a))
plot(real(alpha1b(5:end)),imag(alpha1b(5:end)))
plot(real(alpha1c),imag(alpha1c))
a1=plot(real(alpha1a(end)),imag(alpha1a(end)),'x');
b1=plot(real(alpha1b(end)),imag(alpha1b(end)),'*');
c1=plot(real(alpha1c(end)),imag(alpha1c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Differential Relay, Case #9')

%Problem 3, part c
%Idea: Calculate Ip and Iq for 3PF at m = 0.5 with sweeping RF. Plot alpha
%as function of RF (should be circle) together with restraining zone.
fprintf('\nProblem 3, part c:\n')
Rf = 0:1:1000;
Zx1 = Zs1+ZL1B12+0.5*ZL1B23;
Zy1 = ZR1+0.5*ZL1B23;
Z1 = Zx1*Zy1/(Zx1+Zy1);
C1 = Zy1/(Zx1+Zy1);
Eth = (1-C1)*ER + C1*Es;
If = Eth./(Z1+Rf);
Ipf = (Es-ER)/(Zs1+ZL1B12+ZL1B23+ZR1);
alpha = Ipf./If;
alphad = (alpha+C1)./(alpha+C1-1);
%Determine maximum possible fault resistance that relay can recognize
for n = 1:length(alphad)
    if (abs(alphad(n))<6) && ((angle(alphad(n))*180/pi)>-B)
        break
    end
end
fprintf('Maximum possible fault resistance for relay to operate: %d ohms\n',Rf(n))
fprintf('\nSee attached plot of alpha as a function of fault resistance.')
figure
circ=plot(real(alphad),imag(alphad),'--');
hold on
grid on
xlim([-12,12]);
ylim([-12,12]);
re=plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
a=plot(real(alphad(1)),imag(alphad(1)),'ok');
b=plot(real(alphad(end)),imag(alphad(end)),'xk');
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Alpha as function of fault resistance (3LG fault at m = 0.5)')
c=plot(real(alphad(n)),imag(alphad(n)),'*');
legend([re,circ,a,c,b],{'Restraining Zone','Alpha','Rf = 0 ohms','Rf = 523 ohms','Rf = 1000 ohms'})

%%
% <include>DigiRelay.m</include>
%
%%
% <include>DigiDiffRelay.m</include>
%