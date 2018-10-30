%Leighton Clarke
%Advanced Protections, HW3
%April 24, 2018
%This program reads fault voltages and currents from a datafile, calculates
%the phasor values, and simulates a digital MHO relay with both
%self-polarized and cross-polarized schemes.
close all
clear all
clc

data = csvread('prot2_hw2_2108.csv');
j = sqrt(-1);
theta = 0:0.01:(2*pi);
VTR = 525000/115;
CTR = 800/5;
ZTR = VTR/CTR;
kk = 180/pi;

%%Problem 1
%Calculate phasors using traditional method
E=525000/sqrt(3);
Zs1=2+30*j;
Zs0=Zs1;
zL1=0.073+0.8*j;
zL0=0.1+2.3*j;
I1=E/(Zs1+zL1*20+Zs1+zL1*20+Zs0+zL0*20);
I2=I1;I0=I1;
V1=E-Zs1*I1;
V2=-Zs1*I2;
V0=-Zs0*I0;
a=1*exp(j*120/kk);
Iar=(I0+I1+I2)/CTR;
Ibr=(I0+(a^2)*I1+a*I2)/CTR;
Icr=(I0+a*I1+(a^2)*I2)/CTR;
Var=(V0+V1+V2)/VTR;
Vbr=(V0+(a^2)*V1+a*V2)/VTR;
Vcr=(V0+a*V1+(a^2)*V2)/VTR;
Ires=Iar+Ibr+Icr;

%MHO relay characteristic
ZL1 = (0.073+j*0.8)*100;
Zr = ZL1 * 0.8;
ZL1s = ZL1/ZTR;
Zrs = Zr/ZTR;
r = abs(Zrs)/2;
rR = r*cos(theta) + real(Zrs)/2;
rX = r*sin(theta) + imag(Zrs)/2;

%Calculate impedances
Zab = (Var-Vbr)/(Iar-Ibr);
Zbc = (Vbr-Vcr)/(Ibr-Icr);
Zca = (Vcr-Var)/(Icr-Iar);
k0 = (zL0-zL1)/(3*zL1);
Zag = Var/(Iar+k0*Ires);
Zbg = Vbr/(Ibr+k0*Ires);
Zcg = Vcr/(Icr+k0*Ires);

%Print results
fprintf('Traditionally calculated values:\n')
fprintf('Va: %f V, angle %f deg\n', abs(Var), angle(Var)*kk)
fprintf('Vb: %f V, angle %f deg\n', abs(Vbr), angle(Vbr)*kk)
fprintf('Vc: %f V, angle %f deg\n', abs(Vcr), angle(Vcr)*kk)
fprintf('Ia: %f A, angle %f deg\n', abs(Iar), angle(Iar)*kk)
fprintf('Ib: %f A, angle %f deg\n', abs(Ibr), angle(Ibr)*kk)
fprintf('Ic: %f A, angle %f deg\n', abs(Icr), angle(Icr)*kk)
fprintf('Zab: %f + j%f secondary ohms\n', real(Zab), imag(Zab))
fprintf('Zbc is undefined\n')
fprintf('Zca: %f + j%f secondary ohms\n', real(Zca), imag(Zca))
fprintf('Zag: %f + j%f secondary ohms\n', real(Zag), imag(Zag))
fprintf('Zbg: %f + j%f secondary ohms\n', real(Zbg), imag(Zbg))
fprintf('Zcg: %f + j%f secondary ohms\n', real(Zcg), imag(Zcg))

%Plot faults and protective zone
figure(1)
plot(rR,rX)
grid on
hold on
plot(real(Zab),imag(Zab),'o')
plot(real(Zca),imag(Zca),'*')
plot(real(Zag),imag(Zag),'xk')
plot(real(Zbg),imag(Zbg),'s')
plot(real(Zcg),imag(Zcg),'d')
plot([0,real(ZL1s)],[0,imag(ZL1s)],'--')
legend('Protective Zone','Zab','Zca','Zag','Zbg','Zcg','ZL1')
legend('Location','south')
xlim([-5,5]);
ylim([-5,5]);
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Problem 1: Self-Polarized MHO Characteristic and Fault Impedances')

%%Problem 2
%Calculate voltage and current phasors
%Read from data file
t0 = data(:,1);
Va = data(:,2)./VTR;
Vb = data(:,3)./VTR;
Vc = data(:,4)./VTR;
Ia = data(:,5)./CTR;
Ib = data(:,6)./CTR;
Ic = data(:,7)./CTR;

%Simulate analog RC filter
fc = 600;
RC = 1 / (2*pi*fc);
Vaf = zeros(length(t0),1);
Vbf = zeros(length(t0),1);
Vcf = zeros(length(t0),1);
Iaf = zeros(length(t0),1);
Ibf = zeros(length(t0),1);
Icf = zeros(length(t0),1);

for k = 2:length(t0)
    Vaf(k)=((t0(k)-t0(k-1))*Va(k)+RC*Vaf(k-1))/((t0(k)-t0(k-1))+RC);
    Vbf(k)=((t0(k)-t0(k-1))*Vb(k)+RC*Vbf(k-1))/((t0(k)-t0(k-1))+RC);
    Vcf(k)=((t0(k)-t0(k-1))*Vc(k)+RC*Vcf(k-1))/((t0(k)-t0(k-1))+RC);
    Iaf(k)=((t0(k)-t0(k-1))*Ia(k)+RC*Iaf(k-1))/((t0(k)-t0(k-1))+RC);
    Ibf(k)=((t0(k)-t0(k-1))*Ib(k)+RC*Ibf(k-1))/((t0(k)-t0(k-1))+RC);
    Icf(k)=((t0(k)-t0(k-1))*Ic(k)+RC*Icf(k-1))/((t0(k)-t0(k-1))+RC);
end

%Interpolate data with lower sample rate
Nsc = 16;
fs = Nsc * 60;
Ts = 1 /fs;
ts = 0:Ts:t0(end);
Vafi = interp1(t0, Vaf, ts);
Vbfi = interp1(t0, Vbf, ts);
Vcfi = interp1(t0, Vcf, ts);
Iafi = interp1(t0, Iaf, ts);
Ibfi = interp1(t0, Ibf, ts);
Icfi = interp1(t0, Icf, ts);

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
n = 1:Nsc;
A = (2/Nsc)*cos(2*pi*n/Nsc);
B = (2/Nsc)*sin(2*pi*n/Nsc);
VaRe = filter(A, 1, Vafi);
VbRe = filter(A, 1, Vbfi);
VcRe = filter(A, 1, Vcfi);
IaRe = filter(A, 1, Iafi);
IbRe = filter(A, 1, Ibfi);
IcRe = filter(A, 1, Icfi);
VaIm = filter(B, 1, Vafi);
VbIm = filter(B, 1, Vbfi);
VcIm = filter(B, 1, Vcfi);
IaIm = filter(B, 1, Iafi);
IbIm = filter(B, 1, Ibfi);
IcIm = filter(B, 1, Icfi);

%Calculate phasors
Vap = (1/sqrt(2))*(VaRe+1i*VaIm).*exp(-1i*2*pi*60*ts);
Vbp = (1/sqrt(2))*(VbRe+1i*VbIm).*exp(-1i*2*pi*60*ts);
Vcp = (1/sqrt(2))*(VcRe+1i*VcIm).*exp(-1i*2*pi*60*ts);
Iap = (1/sqrt(2))*(IaRe+1i*IaIm).*exp(-1i*2*pi*60*ts);
Ibp = (1/sqrt(2))*(IbRe+1i*IbIm).*exp(-1i*2*pi*60*ts);
Icp = (1/sqrt(2))*(IcRe+1i*IcIm).*exp(-1i*2*pi*60*ts);

%Create arrays of apparent impedances as seen by relay
Zabr = (Vap-Vbp)./(Iap-Ibp);
Zbcr = (Vbp-Vcp)./(Ibp-Icp);
Zcar = (Vcp-Vap)./(Icp-Iap);
Iresr = Iap+Ibp+Icp;
Zagr = Vap./(Iap+k0*Iresr);
Zbgr = Vbp./(Ibp+k0*Iresr);
Zcgr = Vcp./(Icp+k0*Iresr);

%Print results
fprintf('\nValues measured by digital relay:\n')
fprintf('Va: %f V, angle %f deg\n', abs(Vap(end)), angle(Vap(end))*kk)
fprintf('Vb: %f V, angle %f deg\n', abs(Vbp(end)), angle(Vbp(end))*kk)
fprintf('Vc: %f V, angle %f deg\n', abs(Vcp(end)), angle(Vcp(end))*kk)
fprintf('Ia: %f A, angle %f deg\n', abs(Iap(end)), angle(Iap(end))*kk)
fprintf('Ib: %f A, angle %f deg\n', abs(Ibp(end)), angle(Ibp(end))*kk)
fprintf('Ic: %f A, angle %f deg\n', abs(Icp(end)), angle(Icp(end))*kk)
fprintf('Zab: %f + j%f secondary ohms\n', real(Zabr(end)), imag(Zabr(end)))
fprintf('Zbc: %f + j%f secondary ohms\n', real(Zbcr(end)), imag(Zbcr(end)))
fprintf('Zca: %f + j%f secondary ohms\n', real(Zcar(end)), imag(Zcar(end)))
fprintf('Zag: %f + j%f secondary ohms\n', real(Zagr(end)), imag(Zagr(end)))
fprintf('Zbg: %f + j%f secondary ohms\n', real(Zbgr(end)), imag(Zbgr(end)))
fprintf('Zcg: %f + j%f secondary ohms\n', real(Zcgr(end)), imag(Zcgr(end)))

%Plot evolution of apparent impedance Zag in complex plane
figure(2)
plot(real(Zagr),imag(Zagr))
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Problem 2: Evolution of apparent impedance Zag')

%Plot steady-state of Zag in MHO characteristic
figure(3)
plot(rR,rX)
hold on
plot(real(Zabr(end)),imag(Zabr(end)),'o')
plot(real(Zcar(end)),imag(Zcar(end)),'*')
plot(real(Zagr(end)),imag(Zagr(end)),'xk')
plot(real(Zbgr(end)),imag(Zbgr(end)),'s')
plot(real(Zcgr(end)),imag(Zcgr(end)),'d')
plot([0,real(ZL1s)],[0,imag(ZL1s)],'--')
grid on
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Problem 2: Steady-State MHO Characteristic of Self-Polarized Relay')
legend('Protective Zone','Zab','Zca','Zag','Zbg','Zcg','ZL1')
legend('Location','south')

%Create digital MHO variables for each datapoint
Sab = (Iap-Ibp)*Zrs - (Vap-Vbp);
Sbc = (Ibp-Icp)*Zrs - (Vbp-Vcp);
Sca = (Icp-Iap)*Zrs - (Vcp-Vap);
Sag = (Iap+k0*Iresr)*Zrs - Vap;
Sbg = (Ibp+k0*Iresr)*Zrs - Vbp;
Scg = (Icp+k0*Iresr)*Zrs - Vcp;

MHOAB = zeros(1,length(Sab));
MHOBC = zeros(1,length(Sab));
MHOCA = zeros(1,length(Sab));
MHOAG = zeros(1,length(Sab));
MHOBG = zeros(1,length(Sab));
MHOCG = zeros(1,length(Sab));

for m = 1:length(MHOAB)
    if real(Sab(m)*conj(Vap(m)-Vbp(m))) > 0
        MHOAB(m) = 1;
    end
    if real(Sbc(m)*conj(Vbp(m)-Vcp(m))) > 0
        MHOBC(m) = 1;
    end
    if real(Sca(m)*conj(Vcp(m)-Vap(m))) > 0
        MHOCA(m) = 1;
    end
    if real(Sag(m)*conj(Vap(m))) > 0
        MHOAG(m) = 1;
    end
    if real(Sbg(m)*conj(Vbp(m))) > 0
        MHOBG(m) = 1;
    end
    if real(Scg(m)*conj(Vcp(m))) > 0
        MHOCG(m) = 1;
    end
end

%Plot digital operation flags
figure('position',[200 200 700 700])
subplot(5,2,[1,2])
%ylabel('Original Voltages')
plot(t0,Va)
title({'Comparison of Original Voltages/Currents and Self-Polarized MHO Variables','(1 to operate)'})
hold on
plot(t0,Vb)
plot(t0,Vc)
ylabel('Voltage')
subplot(5,2,[3,4])
plot(t0,Ia)
hold on
plot(t0,Ib)
plot(t0,Ic)
ylabel('Current')
subplot(5,2,5)
plot(ts,MHOAB)
ylabel('MHOAB')
subplot(5,2,6)
plot(ts,MHOBC)
ylabel('MHOBC')
subplot(5,2,7)
plot(ts,MHOCA)
ylabel('MHOCA')
subplot(5,2,8)
plot(ts,MHOAG)
ylabel('MHOAG')
subplot(5,2,9)
plot(ts,MHOBG)
ylabel('MHOBG')
subplot(5,2,10)
plot(ts,MHOCG)
ylabel('MHOCG')

%%Problem 3
%Calculate cross-polarized MHO characteristic for traditionally calculated
%values (AG relay only)
Vpolr = (Vbr-Vcr)*exp(j*90/kk);
Vpr = Var - Vpolr;
Zpr = Vpr/(Iar+k0*Ires);
C = (Zrs + Zpr)/2;
rcross = abs(Zrs-Zpr)/2;
rRcross = rcross*cos(theta)+real(C);
rXcross = rcross*sin(theta)+imag(C);

%Plot both characteristics and apparent impedance Zag
figure(5)
plot(rRcross,rXcross)
hold on
plot(rR,rX,'--')
plot(real(Zag),imag(Zag),'xk')
plot([0,real(Zrs)],[0,imag(Zrs)],':k')
plot([0,real(Zpr)],[0,imag(Zpr)],'-.k')
legend('Cross-polarized characteristic','Self-polarized characteristic','Zag','Zr','Zp')
xlim([-3,3])
ylim([-2,4])
grid on
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Problem 3: Comparison of Self-Polarized and Cross-Polarized MHO Characteristics')

%Calculate MHO operation variables using cross-polarizing scheme
MHOABcross = zeros(1,length(Sab));
MHOBCcross = zeros(1,length(Sab));
MHOCAcross = zeros(1,length(Sab));
MHOAGcross = zeros(1,length(Sab));
MHOBGcross = zeros(1,length(Sab));
MHOCGcross = zeros(1,length(Sab));

for m = 1:length(MHOABcross)
    if real(Sab(m)*conj(Vcp(m)*exp(-j*90/kk))) > 0
        MHOABcross(m) = 1;
    end
    if real(Sbc(m)*conj(Vap*exp(-j*90/kk))) > 0
        MHOBCcross(m) = 1;
    end
    if real(Sca(m)*conj(Vbp(m)*exp(-j*90/kk))) > 0
        MHOCAcross(m) = 1;
    end
    if real(Sag(m)*conj((Vbp(m)-Vcp(m))*exp(j*90/kk))) > 0
        MHOAGcross(m) = 1;
    end
    if real(Sbg(m)*conj((Vcp(m)-Vap(m))*exp(j*90/kk))) > 0
        MHOBGcross(m) = 1;
    end
    if real(Scg(m)*conj((Vap(m)-Vbp(m))*exp(j*90/kk))) > 0
        MHOCGcross(m) = 1;
    end
end

%Plot MHO variables for cross-polarizing scheme
figure('position',[200 200 700 700])
subplot(5,2,[1,2])
plot(t0,Va)
title({'Comparison of Original Voltages/Currents and Cross-Polarized MHO Variables','(1 to operate)'})
hold on
plot(t0,Vb)
plot(t0,Vc)
ylabel('Voltage')
subplot(5,2,[3,4])
plot(t0,Ia)
hold on
plot(t0,Ib)
plot(t0,Ic)
ylabel('Current')
subplot(5,2,5)
plot(ts,MHOABcross)
ylabel('MHOAB')
subplot(5,2,6)
plot(ts,MHOBCcross)
ylabel('MHOBC')
subplot(5,2,7)
plot(ts,MHOCAcross)
ylabel('MHOCA')
subplot(5,2,8)
plot(ts,MHOAGcross)
ylabel('MHOAG')
subplot(5,2,9)
plot(ts,MHOBGcross)
ylabel('MHOBG')
subplot(5,2,10)
plot(ts,MHOCGcross)
ylabel('MHOCG')

%Plot AG element MHO characteristics and final value of Zag at relay
Vpolp = (Vbp-Vcp).*exp(j*90/kk);
Vpp = Vap - Vpolp;
Zpp = Vpp ./ (Iap+k0*Iresr);
Cp = (Zrs + Zpp)/2;
rcrossp = abs(Zrs-Zpp)/2;
rRcrossp = rcrossp(end)*cos(theta)+real(Cp(end));
rXcrossp = rcrossp(end)*sin(theta)+imag(Cp(end));
figure(7)
plot(rRcrossp,rXcrossp)
hold on
plot(rR,rX,'--')
plot(real(Zagr(end)),imag(Zagr(end)),'kx')
grid on
xlim([-3,3])
ylim([-2,4])
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Steady-State Cross-Polarized AG Relay (As Measured)')
legend('Cross-polarized characteristic','Self-polarized characteristic','Zag')
