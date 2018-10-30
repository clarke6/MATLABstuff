%Leighton Clarke
%Advanced Protections, Exam 1
%May 13, 2018

%This program was used to solve problem 2. There are three parts, as
%indicated by the comments.
%%
% <include>DigiRelay2.m</include>
%
%%Part a
%Read data and define system parameters from given information
close all
fault1 = csvread('fault01_PROB2.csv');
fault2 = csvread('fault02_PROB2.csv');
CTR = 800/5;
VTR = 500000/sqrt(3)/67;
ZTR = VTR/CTR;
j = sqrt(-1);
theta = 0:0.01:(2*pi);
kk = 180/pi;
zL0 = 0.1+j*2.6;
zL1 = 0.073+j*0.8;
k0 = (zL0-zL1)/(3*zL1);
ZL1 = zL1*100;
Zr = ZL1*.8/ZTR;
%Run original data through digital relay filtering algorithm
[ts1, va1, vb1, vc1, ia1, ib1, ic1] = DigiRelay2(fault1,CTR,VTR);
[ts2, va2, vb2, vc2, ia2, ib2, ic2] = DigiRelay2(fault2,CTR,VTR);
ires1 = ia1+ib1+ic1;
ires2 = ia2+ib2+ic2;

%Plot original voltage and current signals
figure(1)
subplot(2,1,1);
plot(fault1(:,1),fault1(:,2))
hold on
plot(fault1(:,1),fault1(:,3))
plot(fault1(:,1),fault1(:,4))
legend({'Va','Vb','Vc'})
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Fault 1 Original Voltage Signal')
subplot(2,1,2);
plot(fault1(:,1),fault1(:,5))
hold on
plot(fault1(:,1),fault1(:,6))
plot(fault1(:,1),fault1(:,7))
legend({'Ia','Ib','Ic'})
xlabel('Time (s)')
ylabel('Current (A)')
title('Fault 1 Original Current Signal')

figure(2)
subplot(2,1,1);
plot(fault2(:,1),fault2(:,2))
hold on
plot(fault2(:,1),fault2(:,3))
plot(fault2(:,1),fault2(:,4))
legend({'Va','Vb','Vc'})
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Fault 2 Original Voltage Signal')
subplot(2,1,2);
plot(fault2(:,1),fault2(:,5))
hold on
plot(fault2(:,1),fault2(:,6))
plot(fault2(:,1),fault2(:,7))
legend({'Ia','Ib','Ic'})
xlabel('Time (s)')
ylabel('Current (A)')
title('Fault 2 Original Current Signal')

%Calculate apparent impedances
Zag1 = va1./(ia1+k0*ires1);
Zbg1 = vb1./(ib1+k0*ires1);
Zcg1 = vc1./(ic1+k0*ires1);
Zab1 = (va1-vb1)./(ia1-ib1);
Zbc1 = (vb1-vc1)./(ib1-ic1);
Zca1 = (vc1-va1)./(ic1-ia1);
Zag2 = va2./(ia2+k0*ires2);
Zbg2 = vb2./(ib2+k0*ires2);
Zcg2 = vc2./(ic2+k0*ires2);
Zab2 = (va2-vb2)./(ia2-ib2);
Zbc2 = (vb2-vc2)./(ib2-ic2);
Zca2 = (vc2-va2)./(ic2-ia2);

%Define polarizing signals using specified cross-polarized scheme
Vpolag1 = (vb1-vc1).*exp(j*90/kk);
Vpolbg1 = (vc1-va1).*exp(j*90/kk);
Vpolcg1 = (va1-vb1).*exp(j*90/kk);
Vpolab1 = vc1.*exp(-j*90/kk);
Vpolbc1 = va1.*exp(-j*90/kk);
Vpolca1 = vb1.*exp(-j*90/kk);
Vpolag2 = (vb2-vc2).*exp(j*90/kk);
Vpolbg2 = (vc2-va2).*exp(j*90/kk);
Vpolcg2 = (va2-vb2).*exp(j*90/kk);
Vpolab2 = vc2.*exp(-j*90/kk);
Vpolbc2 = va2.*exp(-j*90/kk);
Vpolca2 = vb2.*exp(-j*90/kk);

%Calculate Vp for each relay
Vpag1 = va1-Vpolag1;
Vpbg1 = vb1-Vpolbg1;
Vpcg1 = vc1-Vpolcg1;
Vpab1 = (va1-vb1)-Vpolab1;
Vpbc1 = (vb1-vc1)-Vpolbc1;
Vpca1 = (vc1-va1)-Vpolca1;
Vpag2 = va2-Vpolag2;
Vpbg2 = vb2-Vpolbg2;
Vpcg2 = vc2-Vpolcg2;
Vpab2 = (va2-vb2)-Vpolab2;
Vpbc2 = (vb2-vc2)-Vpolbc2;
Vpca2 = (vc2-va2)-Vpolca2;

%Calculate Zp for each relay
Zpag1 = Vpag1./(ia1+k0*ires1);
Zpbg1 = Vpbg1./(ib1+k0*ires1);
Zpcg1 = Vpcg1./(ic1+k0*ires1);
Zpab1 = Vpab1./(ia1-ib1);
Zpbc1 = Vpbc1./(ib1-ic1);
Zpca1 = Vpca1./(ic1-ia1);
Zpag2 = Vpag2./(ia2+k0*ires2);
Zpbg2 = Vpbg2./(ib2+k0*ires2);
Zpcg2 = Vpcg2./(ic2+k0*ires2);
Zpab2 = Vpab2./(ia2-ib2);
Zpbc2 = Vpbc2./(ib2-ic2);
Zpca2 = Vpca2./(ic2-ia2);

%Calculate Zc (center impedance) for each relay
Zcag1 = (Zr+Zpag1)/2;
Zcbg1 = (Zr+Zpbg1)/2;
Zccg1 = (Zr+Zpcg1)/2;
Zcab1 = (Zr+Zpab1)/2;
Zcbc1 = (Zr+Zpbc1)/2;
Zcca1 = (Zr+Zpca1)/2;
Zcag2 = (Zr+Zpag2)/2;
Zcbg2 = (Zr+Zpbg2)/2;
Zccg2 = (Zr+Zpcg2)/2;
Zcab2 = (Zr+Zpab2)/2;
Zcbc2 = (Zr+Zpbc2)/2;
Zcca2 = (Zr+Zpca2)/2;

%Calculate radius of each relay characteristic
rag1 = abs(Zr-Zpag1)/2;
rbg1 = abs(Zr-Zpbg1)/2;
rcg1 = abs(Zr-Zpcg1)/2;
rab1 = abs(Zr-Zpab1)/2;
rbc1 = abs(Zr-Zpbc1)/2;
rca1 = abs(Zr-Zpca1)/2;
rag2 = abs(Zr-Zpag2)/2;
rbg2 = abs(Zr-Zpbg2)/2;
rcg2 = abs(Zr-Zpcg2)/2;
rab2 = abs(Zr-Zpab2)/2;
rbc2 = abs(Zr-Zpbc2)/2;
rca2 = abs(Zr-Zpca2)/2;

%Define X and Y components of relay characteristics for plotting
ag1X = rag1(end)*cos(theta)+real(Zcag1(end));
ag1Y = rag1(end)*sin(theta)+imag(Zcag1(end));
bg1X = rbg1(end)*cos(theta)+real(Zcbg1(end));
bg1Y = rbg1(end)*sin(theta)+imag(Zcbg1(end));
cg1X = rcg1(end)*cos(theta)+real(Zccg1(end));
cg1Y = rcg1(end)*sin(theta)+imag(Zccg1(end));
ab1X = rab1(end)*cos(theta)+real(Zcab1(end));
ab1Y = rab1(end)*sin(theta)+imag(Zcab1(end));
bc1X = rbc1(end)*cos(theta)+real(Zcbc1(end));
bc1Y = rbc1(end)*sin(theta)+imag(Zcbc1(end));
ca1X = rca1(end)*cos(theta)+real(Zcca1(end));
ca1Y = rca1(end)*sin(theta)+imag(Zcca1(end));

%Plot each relay's final characteristic together with evolution of apparent
%impedances over time
figure(3)
plot(ag1X,ag1Y);
hold on
grid on
plot(real(Zag1),imag(Zag1));
title('Zag for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpag1(end))],[0,imag(Zpag1(end))],'--');
legend([a,b],{'Zr','Zp'});
figure(4)
plot(bg1X,bg1Y);
hold on
grid on
plot(real(Zbg1),imag(Zbg1));
title('Zbg for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpbg1(end))],[0,imag(Zpbg1(end))],'--');
legend([a,b],{'Zr','Zp'});
figure(5)
plot(cg1X,cg1Y);
hold on
grid on
plot(real(Zcg1),imag(Zcg1));
title('Zcg for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpcg1(end))],[0,imag(Zpcg1(end))],'--');
legend([a,b],{'Zr','Zp'});
figure(6)
plot(ab1X,ab1Y);
hold on
grid on
plot(real(Zab1),imag(Zab1));
title('Zab for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpab1(end))],[0,imag(Zpab1(end))],'--');
legend([a,b],{'Zr','Zp'});
figure(7)
plot(bc1X,bc1Y);
hold on
grid on
plot(real(Zbc1),imag(Zbc1));
title('Zbc for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpbc1(end))],[0,imag(Zpbc1(end))],'--');
legend([a,b],{'Zr','Zp'});
figure(8)
plot(ca1X,ca1Y);
hold on
grid on
plot(real(Zca1),imag(Zca1));
title('Zca for Fault 1')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)]);
b=plot([0,real(Zpca1(end))],[0,imag(Zpca1(end))],'--');
legend([a,b],{'Zr','Zp'});

%Repeat for fault 2
ag2X = rag2(end)*cos(theta)+real(Zcag2(end));
ag2Y = rag2(end)*sin(theta)+imag(Zcag2(end));
bg2X = rbg2(end)*cos(theta)+real(Zcbg2(end));
bg2Y = rbg2(end)*sin(theta)+imag(Zcbg2(end));
cg2X = rcg2(end)*cos(theta)+real(Zccg2(end));
cg2Y = rcg2(end)*sin(theta)+imag(Zccg2(end));
ab2X = rab2(end)*cos(theta)+real(Zcab2(end));
ab2Y = rab2(end)*sin(theta)+imag(Zcab2(end));
bc2X = rbc2(end)*cos(theta)+real(Zcbc2(end));
bc2Y = rbc2(end)*sin(theta)+imag(Zcbc2(end));
ca2X = rca2(end)*cos(theta)+real(Zcca2(end));
ca2Y = rca2(end)*sin(theta)+imag(Zcca2(end));


figure(9)
plot(ag2X,ag2Y)
hold on
grid on
plot(real(Zag2),imag(Zag2));
title('Zag for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpag2(end))],[0,imag(Zpag2(end))],'--k');
legend([a,b],{'Zr','Zp'});
figure(10)
plot(bg2X,bg2Y);
hold on
grid on
plot(real(Zbg2),imag(Zbg2));
title('Zbg for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpbg2(end))],[0,imag(Zpbg2(end))],'--k');
legend([a,b],{'Zr','Zp'});
figure(11)
plot(cg2X,cg2Y);
hold on
grid on
plot(real(Zcg2),imag(Zcg2));
title('Zcg for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpcg2(end))],[0,imag(Zpcg2(end))],'--k');
legend([a,b],{'Zr','Zp'});
figure(12)
plot(ab2X,ab2Y);
hold on
grid on
plot(real(Zab2),imag(Zab2));
title('Zab for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpab2(end))],[0,imag(Zpab2(end))],'--k');
legend([a,b],{'Zr','Zp'});
figure(13)
plot(bc2X,bc2Y);
hold on
grid on
plot(real(Zbc2),imag(Zbc2));
title('Zbc for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpbc2(end))],[0,imag(Zpbc2(end))],'--k');
legend([a,b],{'Zr','Zp'});
figure(14)
plot(ca2X,ca2Y);
hold on
grid on
plot(real(Zca2),imag(Zca2));
title('Zca for Fault 2')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)') 
a=plot([0,real(Zr)],[0,imag(Zr)],'k');
b=plot([0,real(Zpca2(end))],[0,imag(Zpca2(end))],'--k');
legend([a,b],{'Zr','Zp'});

%Create digital MHO variables for each datapoint (Fault 1)
S1ab1 = (ia1-ib1)*Zr - (va1-vb1);
S1bc1 = (ib1-ic1)*Zr - (vb1-vc1);
S1ca1 = (ic1-ia1)*Zr - (vc1-va1);
S1ag1 = (ia1+k0*ires1)*Zr - va1;
S1bg1 = (ib1+k0*ires1)*Zr - vb1;
S1cg1 = (ic1+k0*ires1)*Zr - vc1;

MHOAB1 = zeros(1,length(ts1));
MHOBC1 = zeros(1,length(ts1));
MHOCA1 = zeros(1,length(ts1));
MHOAG1 = zeros(1,length(ts1));
MHOBG1 = zeros(1,length(ts1));
MHOCG1 = zeros(1,length(ts1));

for m = 1:length(ts1)
    if real(S1ab1(m)*conj(Vpolab1(m))) > 0
        MHOAB1(m) = 1;
    end
    if real(S1bc1(m)*conj(Vpolbc1(m))) > 0
        MHOBC1(m) = 1;
    end
    if real(S1ca1(m)*conj(Vpolca1(m))) > 0
        MHOCA1(m) = 1;
    end
    if real(S1ag1(m)*conj(Vpolag1(m))) > 0
        MHOAG1(m) = 1;
    end
    if real(S1bg1(m)*conj(Vpolbg1(m))) > 0
        MHOBG1(m) = 1;
    end
    if real(S1cg1(m)*conj(Vpolcg1(m))) > 0
        MHOCG1(m) = 1;
    end
end

figure(15)
subplot(3,2,1)
plot(ts1,MHOAB1);
ylabel('MHOAB')
ylim([-.25 1.25])
subplot(3,2,2)
plot(ts1,MHOBC1);
ylabel('MHOBC')
ylim([-.25 1.25])
subplot(3,2,3)
plot(ts1,MHOCA1);
ylabel('MHOCA')
ylim([-.25 1.25])
subplot(3,2,4)
plot(ts1,MHOAG1);
ylabel('MHOAG')
ylim([-.25 1.25])
subplot(3,2,5)
plot(ts1,MHOBG1);
ylabel('MHOBG')
xlabel('Time (s)')
ylim([-.25 1.25])
subplot(3,2,6)
plot(ts1,MHOCG1);
ylabel('MHOCG')
xlabel('Time (s)')
ylim([-.25 1.25])
suptitle('MHO Digital Outputs (Fault 1)')

%Create digital MHO variables for each datapoint (Fault 2)
S1ab2 = (ia2-ib2)*Zr - (va2-vb2);
S1bc2 = (ib2-ic2)*Zr - (vb2-vc2);
S1ca2 = (ic2-ia2)*Zr - (vc2-va2);
S1ag2 = (ia2+k0*ires2)*Zr - va2;
S1bg2 = (ib2+k0*ires2)*Zr - vb2;
S1cg2 = (ic2+k0*ires2)*Zr - vc2;

MHOAB2 = zeros(1,length(ts2));
MHOBC2 = zeros(1,length(ts2));
MHOCA2 = zeros(1,length(ts2));
MHOAG2 = zeros(1,length(ts2));
MHOBG2 = zeros(1,length(ts2));
MHOCG2 = zeros(1,length(ts2));

for m = 1:length(ts2)
    if real(S1ab2(m)*conj(Vpolab2(m))) > 0
        MHOAB2(m) = 1;
    end
    if real(S1bc2(m)*conj(Vpolbc2(m))) > 0
        MHOBC2(m) = 1;
    end
    if real(S1ca2(m)*conj(Vpolca2(m))) > 0
        MHOCA2(m) = 1;
    end
    if real(S1ag2(m)*conj(Vpolag2(m))) > 0
        MHOAG2(m) = 1;
    end
    if real(S1bg2(m)*conj(Vpolbg2(m))) > 0
        MHOBG2(m) = 1;
    end
    if real(S1cg2(m)*conj(Vpolcg2(m))) > 0
        MHOCG2(m) = 1;
    end
end

figure(16)
subplot(3,2,1)
plot(ts2,MHOAB2);
ylabel('MHOAB')
ylim([-.25 1.25])
subplot(3,2,2)
plot(ts2,MHOBC2);
ylabel('MHOBC')
ylim([-.25 1.25])
subplot(3,2,3)
plot(ts2,MHOCA2);
ylabel('MHOCA')
ylim([-.25 1.25])
subplot(3,2,4)
plot(ts2,MHOAG2);
ylabel('MHOAG')
ylim([-.25 1.25])
subplot(3,2,5)
plot(ts2,MHOBG2);
ylabel('MHOBG')
xlabel('Time (s)')
ylim([-.25 1.25])
subplot(3,2,6)
plot(ts2,MHOCG2);
ylabel('MHOCG')
ylim([-.25 1.25])
xlabel('Time (s)')
suptitle('MHO Digital Outputs (Fault 2)')

%%Part b
%Calculate and plot m using modified Takagi method (let kc = 1)
%Fault 1 is b-g, fault 2 is c-g.
m1 = imag(Zbg1)*ZTR./imag(ZL1);
figure
plot(ts1,m1)
title('m vs Time for Fault 1')
xlabel('Time (s)')
ylabel('m')
dist1 = m1(end)*100;
fprintf('Fault 1 distance: %.4g miles\n', dist1)

m2 = imag(Zcg2)*ZTR./imag(ZL1);
figure
plot(ts2,m2)
title('m vs Time for Fault 2')
xlabel('Time (s)')
ylabel('m')
dist2 = m2(end)*100;
fprintf('Fault 2 distance: %.4g miles\n', dist2)

%%Part c
%calculate total inductance and capacitance of lines B1-B2 and B2-B3
XB1B2 = 50*0.8;
LB1B2 = XB1B2/(2*pi*60);
LB1B2 = LB1B2 / (50*1609.34); %Convert to per meter quantities
CB1B2 = 0.013e-6*50 / (50*1609.34);
XB2B3 = 100*0.8;
LB2B3 = XB2B3/(2*pi*60);
LB2B3 = LB2B3 / (100*1609.34); %Convert to per meter quantities;
CB2B3 = 0.013e-6*100 / (100*1609.34);
%Calculate velocity of traveling wave on each line (should be equal)
vB1B2 = 1/sqrt(LB1B2*CB1B2);
vB2B3 = 1/sqrt(LB2B3*CB2B3);
%Calculate time difference for each fault
dt1 = (m1(end)*100*1609.34*2 - 100*1609.34)/vB2B3;
dt2 = ((50+m2(end)*100)*1609.34*2 - 50*1609.34)/vB1B2;
fprintf('Time difference between traveling wave arrival at buses B2 and B3 for fault 1: %f s\n', dt1)
fprintf('Time difference between traveling wave arrival at buses B1 and B2 for fault 2: %f s\n', dt2)

