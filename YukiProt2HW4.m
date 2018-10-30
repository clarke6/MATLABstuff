%Yuki Kato
%EE537 HW4-5

%%
% <include>RelayFilters.m</include>
%
%%
% <include>DifRelFil.m</include>
%
%%
% <include>agfault.m</include>
%
%%
% <include>bcfault.m</include>
%
%%
% <include>bcgfault.m</include>
%

CTR=800/5;
VTR=500000/sqrt(3)/67;
j=sqrt(-1);
kk = 180/pi;
Zs1 = 1 + j*10;
Zs0 = 2 + j*30;
Zr1 = j*20;
Zr0 = j*10;
zL1 = 0.073 + j*0.8;
zL0 = 0.1 + j*2.6;
ZL1_23 = 100*zL1;
ZL1_12 = 50*zL1;
ZL0_23 = zL0*100;
ZL0_12 = zL0*50;
Zs1F1 = ZL1_12+Zs1;
Zs0F1 = ZL0_12+Zs0;
ZL1F1 = ZL1_23;
ZL0F1 = ZL0_23;
Zr1F2 = ZL1_23+Zr1;
Zr0F2 = ZL0_23+Zr0;
ZL1F2 = ZL1_12;
ZL0F2 = ZL0_12;
V = 500000;
Es = 1.002*exp(j*12.5/kk)*V/sqrt(3);
Er = 0.994*V/sqrt(3);

fprintf('\n------Problem 1, part a------\n')
fprintf('\n---Phasors Calculated Using Sequence Components---\n')
fprintf('\nCase 1:\n')
agfault(Es,Zs1F1,Zs0F1,Er,Zr1,Zr0,ZL1F1,ZL0F1,0.2,0,CTR,VTR)
fprintf('\nCase 2:\n')
agfault(Es,Zs1F1,Zs0F1,Er,Zr1,Zr0,ZL1F1,ZL0F1,0.2,2,CTR,VTR)
fprintf('\nCase 3:\n')
bcfault(Es,Zs1F1,Zs0F1,Er,Zr1,Zr0,ZL1F1,ZL0F1,0.2,0,CTR,VTR)
fprintf('\nCase 4:\n')
bcfault(Es,Zs1F1,Zs0F1,Er,Zr1,Zr0,ZL1F1,ZL0F1,0.2,2,CTR,VTR)
fprintf('\nCase 5:\n')
agfault(Er,Zr1F2,Zr0F2,Es,Zs1,Zs0,ZL1F2,ZL0F2,0.2,0,CTR,VTR)
fprintf('\nCase 6:\n')
agfault(Er,Zr1F2,Zr0F2,Es,Zs1,Zs0,ZL1F2,ZL0F2,0.2,2,CTR,VTR)
fprintf('\nCase 7:\n')
bcfault(Er,Zr1F2,Zr0F2,Es,Zs1,Zs0,ZL1F2,ZL0F2,0.2,0,CTR,VTR)
fprintf('\nCase 8:\n')
bcfault(Er,Zr1F2,Zr0F2,Es,Zs1,Zs0,ZL1F2,ZL0F2,0.2,2,CTR,VTR)
fprintf('\nCase 9:\n')
bcgfault(Es,Zs1F1,Zs0F1,Er,Zr1,Zr0,ZL1F1,ZL0F1,0.2,CTR,VTR)

case1=csvread('f1_ag_r0.csv');
case2=csvread('f1_ag_r2.csv');
case3=csvread('f1_bc_r0.csv');
case4=csvread('f1_bc_r2.csv');
case5=csvread('f2_ag_r0.csv');
case6=csvread('f2_ag_r2.csv');
case7=csvread('f2_bc_r0.csv');
case8=csvread('f2_bc_r2.csv');
case9=csvread('f1_bcg_r0.csv');

fprintf('\n---Phasors from Digital Relay Simulation---\n')
fprintf('\nCase 1:\n')
[Va,Vb,Vc,Ia,Ib,Ic] = RelayFilters(case1,CTR,VTR);
fprintf('\nCase 2:\n')
RelayFilters(case2,CTR,VTR);
fprintf('\nCase 3:\n')
RelayFilters(case3,CTR,VTR);
fprintf('\nCase 4:\n')
RelayFilters(case4,CTR,VTR);
fprintf('\nCase 5:\n')
RelayFilters(case5,CTR,VTR);
fprintf('\nCase 6:\n')
RelayFilters(case6,CTR,VTR);
fprintf('\nCase 7:\n')
RelayFilters(case7,CTR,VTR);
fprintf('\nCase 8:\n')
RelayFilters(case8,CTR,VTR);
fprintf('\nCase 9:\n')
RelayFilters(case9,CTR,VTR);

fprintf('\n------Problem 1, part b------\n')
Ipf = (Es-Er)/(Zs1+ZL1_12+ZL1_23+Zr1);
VB2pf = Es - Ipf*(Zs1F1);
SB2pf = 3*VB2pf*conj(Ipf);
fprintf('Prefault current: %f A, angle %f degrees\n', abs(Ipf), angle(Ipf)*kk)
fprintf('Prefault power exported  at bus 2: %f MW + %f MVAR\n', real(SB2pf)/10^6,imag(SB2pf)/10^6)

fprintf('\n------Problem 1, part c------\n')
fprintf('See attached text.\n')

fprintf('\n------Problem 1, part d------\n')
k0 = (zL0-zL1)/(3*zL1);
ZTR = VTR/CTR;
Ires = Ia+Ib+Ic;
Zag = Va./(Ia+k0*Ires)*ZTR;
mF1 = imag(Zag)./imag(ZL1_23);
fprintf('Calculated Fault Distance: %f miles\n',mF1(end)*100)
fprintf('See plot of m vs time.\n')
t0 = case1(:,1);
t = linspace(0,t0(end),length(mF1));
figure(1)
plot(t,mF1)
xlabel('Time (s)')
ylabel('Per unit distance (m)')
title('Per unit distance vs time')

fprintf('\n------Problem 2------\n')
fprintf('See attached diagram for relay labels\n')
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
fprintf('Relay A:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrAB1),angle(ZrAB1)*kk,abs(ZrAB2),angle(ZrAB2)*kk)
fprintf('Relay B:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrAB1),angle(ZrAB1)*kk,abs(ZrAB2),angle(ZrAB2)*kk)
fprintf('Relay C:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrCD1),angle(ZrCD1)*kk,abs(ZrCD2),angle(ZrCD2)*kk)
fprintf('Relay D:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_12,imag(k0_12),abs(ZrCD1),angle(ZrCD1)*kk,abs(ZrCD2),angle(ZrCD2)*kk)
fprintf('Relay E:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_13,imag(k0_13),abs(ZrE1),angle(ZrE1)*kk,abs(ZrE2),angle(ZrE2)*kk)
fprintf('Relay F:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_14,imag(k0_14),abs(ZrF1),angle(ZrF1)*kk,abs(ZrF2),angle(ZrF2)*kk)
fprintf('Relay G:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_13,imag(k0_13),abs(ZrG1),angle(ZrG1)*kk,abs(ZrG2),angle(ZrG2)*kk)
fprintf('Relay H:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_34,imag(k0_34),abs(ZrH1),angle(ZrH1)*kk,abs(ZrH2),angle(ZrH2)*kk)
fprintf('Relay I:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_34,imag(k0_34),abs(ZrI1),angle(ZrI1)*kk,abs(ZrI2),angle(ZrI2)*kk)
fprintf('Relay J:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_14,imag(k0_14),abs(ZrJ1),angle(ZrJ1)*kk,abs(ZrJ2),angle(ZrJ2)*kk)
fprintf('Relay K:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_45,imag(k0_45),abs(ZrK1),angle(ZrK1)*kk,abs(ZrK2),angle(ZrK2)*kk)
fprintf('Relay L:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_45,imag(k0_45),abs(ZrL1),angle(ZrL1)*kk,abs(ZrL2),angle(ZrL2)*kk)
fprintf('Relay M:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_25,imag(k0_25),abs(ZrM1),angle(ZrM1)*kk,abs(ZrM2),angle(ZrM2)*kk)
fprintf('Relay N:\nk0 = %f+j%f, Zone 1 = %f ohms (angle %f degrees) Zone 2 = %f ohms (angle %f degrees)\n', k0_25,imag(k0_25),abs(ZrN1),angle(ZrN1)*kk,abs(ZrN2),angle(ZrN2)*kk)
fprintf('\nFor all relays, the time delay for zone 2 is 0.3 seconds, and for zone 1 is instantaneous.\n')

fprintf('\n------Problem 3, part a------\n')
case1=csvread('f1agxx.csv');
case2=csvread('f1r2agxx.csv');
case3=csvread('f1bcxx.csv');
case4=csvread('f1r2bcxx.csv');
case5=csvread('f2agxx.csv');
case6=csvread('f2r2agxx.csv');
case7=csvread('f2bcx.csv');
case8=csvread('f2r2bcxx.csv');
case9=csvread('f1bcgxx.csv');
CTR=1200/5;
fprintf('Sending and receiving currents from digital relay simulation:\n')
fprintf('\nCase 1:\n')
[ISA1, ISB1, ISC1, IRA1, IRB1, IRC1] = DifRelFil(case1,CTR);
fprintf('\nCase 2:\n')
[ISA2, ISB2, ISC2, IRA2, IRB2, IRC2] = DifRelFil(case2,CTR);
fprintf('\nCase 3:\n')
[ISA3, ISB3, ISC3, IRA3, IRB3, IRC3] = DifRelFil(case3,CTR);
fprintf('\nCase 4:\n')
[ISA4, ISB4, ISC4, IRA4, IRB4, IRC4] = DifRelFil(case4,CTR);
fprintf('\nCase 5:\n')
[ISA5, ISB5, ISC5, IRA5, IRB5, IRC5] = DifRelFil(case5,CTR);
fprintf('\nCase 6:\n')
[ISA6, ISB6, ISC6, IRA6, IRB6, IRC6] = DifRelFil(case6,CTR);
fprintf('\nCase 7:\n')
[ISA7, ISB7, ISC7, IRA7, IRB7, IRC7] = DifRelFil(case7,CTR);
fprintf('\nCase 8:\n')
[ISA8, ISB8, ISC8, IRA8, IRB8, IRC8] = DifRelFil(case8,CTR);
fprintf('\nCase 9:\n')
[ISA9, ISB9, ISC9, IRA9, IRB9, IRC9] = DifRelFil(case9,CTR);

fprintf('\n------Problem 3, part b------\n')
fprintf('See attached plots for each case.\n')
fprintf('Case 1: "A" phase relay operates.\n')
fprintf('Case 2: "A" phase relay operates.\n')
fprintf('Case 3: "B" and "C" phase relays operate.\n')
fprintf('Case 4: "B" and "C" phase relays operate.\n')
fprintf('Case 5: No relays operate.\n')
fprintf('Case 6: No relays operate.\n')
fprintf('Case 7: No relays operate.\n')
fprintf('Case 8: No relays operate.\n')
fprintf('Case 9: "B" and "C" phase relays operate.\n')

alpha1a = ISA1./IRA1;
alpha1b = ISB1./IRB1;
alpha1c = ISC1./IRC1;
alpha2a = ISA2./IRA2;
alpha2b = ISB2./IRB2;
alpha2c = ISC2./IRC2;
alpha3a = ISA3./IRA3;
alpha3b = ISB3./IRB3;
alpha3c = ISC3./IRC3;
alpha4a = ISA4./IRA4;
alpha4b = ISB4./IRB4;
alpha4c = ISC4./IRC4;
alpha5a = ISA5./IRA5;
alpha5b = ISB5./IRB5;
alpha5c = ISC5./IRC5;
alpha6a = ISA6./IRA6;
alpha6b = ISB6./IRB6;
alpha6c = ISC6./IRC6;
alpha7a = ISA7./IRA7;
alpha7b = ISB7./IRB7;
alpha7c = ISC7./IRC7;
alpha8a = ISA8./IRA8;
alpha8b = ISB8./IRB8;
alpha8c = ISC8./IRC8;
alpha9a = ISA9./IRA9;
alpha9b = ISB9./IRB9;
alpha9c = ISC9./IRC9;

r1 = 1/6;
r2 = 6;
B = 98;
Rtheta = (-B/kk):0.01:(B/kk);

figure(2)
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
xlim([-8 8]);
ylim([-8 8]);
title('Case #1')

figure(3)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha2a),imag(alpha2a))
plot(real(alpha2b(5:end)),imag(alpha2b(5:end)))
plot(real(alpha2c),imag(alpha2c))
a1=plot(real(alpha2a(end)),imag(alpha2a(end)),'x');
b1=plot(real(alpha2b(end)),imag(alpha2b(end)),'*');
c1=plot(real(alpha2c(end)),imag(alpha2c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #2')

figure(4)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha3a),imag(alpha3a))
plot(real(alpha3b(5:end)),imag(alpha3b(5:end)))
plot(real(alpha3c),imag(alpha3c))
a1=plot(real(alpha3a(end)),imag(alpha3a(end)),'x');
b1=plot(real(alpha3b(end)),imag(alpha3b(end)),'*');
c1=plot(real(alpha3c(end)),imag(alpha3c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #3')

figure(5)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha4a),imag(alpha4a))
plot(real(alpha4b(5:end)),imag(alpha4b(5:end)))
plot(real(alpha4c),imag(alpha4c))
a1=plot(real(alpha4a(end)),imag(alpha4a(end)),'x');
b1=plot(real(alpha4b(end)),imag(alpha4b(end)),'*');
c1=plot(real(alpha4c(end)),imag(alpha4c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #4')

figure(6)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha5a),imag(alpha5a))
plot(real(alpha5b(5:end)),imag(alpha5b(5:end)))
plot(real(alpha5c),imag(alpha5c))
a1=plot(real(alpha5a(end)),imag(alpha5a(end)),'x');
b1=plot(real(alpha5b(end)),imag(alpha5b(end)),'*');
c1=plot(real(alpha5c(end)),imag(alpha5c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #5')

figure(7)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha6a),imag(alpha6a))
plot(real(alpha6b(5:end)),imag(alpha6b(5:end)))
plot(real(alpha6c),imag(alpha6c))
a1=plot(real(alpha6a(end)),imag(alpha6a(end)),'x');
b1=plot(real(alpha6b(end)),imag(alpha6b(end)),'*');
c1=plot(real(alpha6c(end)),imag(alpha6c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #6')

figure(8)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha7a),imag(alpha7a))
plot(real(alpha7b(5:end)),imag(alpha7b(5:end)))
plot(real(alpha7c),imag(alpha7c))
a1=plot(real(alpha7a(end)),imag(alpha7a(end)),'x');
b1=plot(real(alpha7b(end)),imag(alpha7b(end)),'*');
c1=plot(real(alpha7c(end)),imag(alpha7c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #7')

figure(9)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha8a),imag(alpha8a))
plot(real(alpha8b(5:end)),imag(alpha8b(5:end)))
plot(real(alpha8c),imag(alpha8c))
a1=plot(real(alpha8a(end)),imag(alpha8a(end)),'x');
b1=plot(real(alpha8b(end)),imag(alpha8b(end)),'*');
c1=plot(real(alpha8c(end)),imag(alpha8c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #8')

figure(10)
plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
hold on
grid on
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
plot(real(alpha9a),imag(alpha9a))
plot(real(alpha9b(5:end)),imag(alpha9b(5:end)))
plot(real(alpha9c),imag(alpha9c))
a1=plot(real(alpha9a(end)),imag(alpha9a(end)),'x');
b1=plot(real(alpha9b(end)),imag(alpha9b(end)),'*');
c1=plot(real(alpha9c(end)),imag(alpha9c(end)),'o');
legend([a1,b1,c1],{'A phase','B phase','C phase'})
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
xlim([-8 8]);
ylim([-8 8]);
title('Case #9')

fprintf('\n------Problem 3, part c------\n')
Rf = 0:1:1000;
Zx1=Zs1+ZL1_12+.5*ZL1_23;
Zy1 = Zr1+0.5*ZL1_23;
Z1 = Zx1*Zy1/(Zx1+Zy1);
C1 = Zy1/(Zx1+Zy1);
Eth = (1-C1)*Er + C1*Es;
If = Eth./(Z1+Rf);
Ipf = (Es-Er)/(Zs1+ZL1_12+ZL1_23+Zr1);
alpha = Ipf./If;
alphad = (alpha+C1)./(alpha+C1-1);
for n = 1:length(alphad)
    if (abs(alphad(n))<6) && ((angle(alphad(n))*180/pi)>-B)
        break
    end
end
fprintf('Maximum possible fault resistance for relay to operate: %d ohms\n',Rf(n))
fprintf('\nSee attached plot of alpha as a function of fault resistance.\n')
figure(11)
circ=plot(real(alphad),imag(alphad),'*');
hold on
grid on
xlim([-12,12]);
ylim([-12,12]);
re=plot(r1*cos(Rtheta),r1*sin(Rtheta),'k');
plot(r2*cos(Rtheta),r2*sin(Rtheta),'k');
plot([r1*cos(B/kk),r2*cos(B/kk)],[r1*sin(B/kk),r2*sin(B/kk)],'k');
plot([r1*cos(-B/kk),r2*cos(-B/kk)],[r1*sin(-B/kk),r2*sin(-B/kk)],'k');
%plot(real(alphad(1)),imag(alphad(1)),'ok');
%plot(real(alphad(end)),imag(alphad(end)),'xk');
xlabel('Re(Ip/Iq)')
ylabel('Im(Ip/Iq)')
title('Effect of Fault Resistance on Differential Relay (3-Phase Fault at location F1)')
%plot(real(alphad(n)),imag(alphad(n)),'*');
legend([re,circ],{'Restraining Zone','Alpha (variable Rf)'})