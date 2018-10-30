%Yuki Kato
%Protections 2, Homework 2
%4-18-2018


clear all
close all


%Read data and create arrays for secondary signals
data = csvread('prot2_hw2_2108.csv');
PTR = 525000/115;
CTR = 800/5;
kk=180/pi;
t0 = data(:,1);
va = data(:,2)./PTR;
vb = data(:,3)./PTR;
vc = data(:,4)./PTR;
ia = data(:,5)./CTR;
ib = data(:,6)./CTR;
ic = data(:,7)./CTR;

%Calculate RC and preallocate memory for new arrays
fc=600;
fc = 600;
RC = 1 / (2*pi*fc);
vaf = zeros(length(t0),1);
vbf = zeros(length(t0),1);
vcf = zeros(length(t0),1);
iaf = zeros(length(t0),1);
ibf = zeros(length(t0),1);
icf = zeros(length(t0),1);
%Use recursive equation to simulate analog low-pass filter
for k = 2:length(t0)
    vaf(k)=((t0(k)-t0(k-1))*va(k)+RC*vaf(k-1))/((t0(k)-t0(k-1))+RC);
    vbf(k)=((t0(k)-t0(k-1))*vb(k)+RC*vbf(k-1))/((t0(k)-t0(k-1))+RC);
    vcf(k)=((t0(k)-t0(k-1))*vc(k)+RC*vcf(k-1))/((t0(k)-t0(k-1))+RC);
    iaf(k)=((t0(k)-t0(k-1))*ia(k)+RC*iaf(k-1))/((t0(k)-t0(k-1))+RC);
    ibf(k)=((t0(k)-t0(k-1))*ib(k)+RC*ibf(k-1))/((t0(k)-t0(k-1))+RC);
    icf(k)=((t0(k)-t0(k-1))*ic(k)+RC*icf(k-1))/((t0(k)-t0(k-1))+RC);
end

%Calculate sampling period and interpolate data
Nsc=16;
fs=Nsc*60;
Dt=1/fs;
t=[0:Dt:t0(length(t0))];
vafr=interp1(t0,vaf,t);
vbfr=interp1(t0,vbf,t);
vcfr=interp1(t0,vcf,t);
iafr=interp1(t0,iaf,t);
ibfr=interp1(t0,ibf,t);
icfr=interp1(t0,icf,t);

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
k=[1:Nsc];
bc=(2/Nsc)*cos(2*pi*k/Nsc);
bs=(2/Nsc)*sin(2*pi*k/Nsc);
vax=filter(bc,1,vafr);
vbx=filter(bc,1,vbfr);
vcx=filter(bc,1,vcfr);
vay=filter(bs,1,vafr);
vby=filter(bs,1,vbfr);
vcy=filter(bs,1,vcfr);
iax=filter(bc,1,iafr);
ibx=filter(bc,1,ibfr);
icx=filter(bc,1,icfr);
iay=filter(bs,1,iafr);
iby=filter(bs,1,ibfr);
icy=filter(bs,1,icfr);

%Create phasor array for each signal
j=sqrt(-1);
Va=(1/sqrt(2))*(vax+j*vay).*exp(-j*2*pi*60*t);
Vb=(1/sqrt(2))*(vbx+j*vby).*exp(-j*2*pi*60*t);
Vc=(1/sqrt(2))*(vcx+j*vcy).*exp(-j*2*pi*60*t);
Ia=(1/sqrt(2))*(iax+j*iay).*exp(-j*2*pi*60*t);
Ib=(1/sqrt(2))*(ibx+j*iby).*exp(-j*2*pi*60*t);
Ic=(1/sqrt(2))*(icx+j*icy).*exp(-j*2*pi*60*t);

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
Var=(V0+V1+V2)/PTR;
Vbr=(V0+(a^2)*V1+a*V2)/PTR;
Vcr=(V0+a*V1+(a^2)*V2)/PTR;

%Create table comparing measured and calculated results
Signal = {'Va'; 'Vb'; 'Vc'; 'Ia'; 'Ib'; 'Ic'};
Relay_V = [abs(Va(end));abs(Vb(end));abs(Vc(end)); ...
    abs(Ia(end));abs(Ib(end));abs(Ic(end))];
Relay_Angle = [0;angle(Vb(end))*kk - angle(Va(end))*kk; ...
    angle(Vc(end))*kk - angle(Va(end))*kk;angle(Ia(end))*kk - angle(Va(end))*kk; ...
    angle(Ib(end))*kk - angle(Va(end))*kk;angle(Ic(end))*kk - angle(Va(end))*kk];
Calculated_V = [abs(Var);abs(Vbr);abs(Vcr);abs(Iar);abs(Ibr);abs(Icr)];
Calculated_Angle = [angle(Var)*kk;angle(Vbr)*kk;angle(Vcr)*kk; ...
    angle(Iar)*kk;angle(Ibr)*kk;angle(Icr)*kk];
T = table(Signal, Relay_V, Relay_Angle, Calculated_V, ...
    Calculated_Angle);
disp(T)



%Create plots
%Phasor angles are plotted both unaltered and then referenced to the angle
%of the a-phase voltage.
figure(1)
subplot(211)
plot(60*t0,va,60*t0,vb,60*t0,vc)
grid on
title('Original instantaneous voltages')
ylabel('Voltage in volts')
xlabel('Time in cycles')
xlim([0 6])

subplot(212)
plot(60*t0,ia,60*t0,ib,60*t0,ic)
grid on
title('Original instantaneous currents')
ylabel('Current in amps')
xlabel('Time in cycles')
xlim([0 6])
set(gca,'YMinorGrid','on')

figure(2);plot(t0,va),grid;hold on;plot(t,abs(Va),'r')
title('Original instantaneous phase-a voltage and magnitude of the phasor')
ylabel('Voltage in volts')
xlabel('Time in seconds')

figure(3);plot(t0,ia),grid;hold on;plot(t,abs(Ia),'r')
title('Original instantaneous phase-a current and magnitude of the phasor')
ylabel('Current in amps')
xlabel('Time in seconds')

figure(4);plot(t,angle(Va)*kk,t,angle(Ia)*kk)
title('Angles in degrees: Blue-->Voltage, Green-->Current')
ylabel('Phase-a voltage and current angles in degrees')
xlabel('Time in seconds')

figure(5)
subplot(3,2,1)
plot(t0,va),grid
ylabel('va')
xlabel('t')
xlim([0 0.1])

subplot(3,2,3)
plot(t0,vaf),grid
ylabel('va after LPF')
xlabel('t')
xlim([0 0.1])

subplot(3,2,5)
plot(t,vax),grid
ylabel('va after COSF')
xlabel('t')
xlim([0 0.1])

subplot(3,2,2)
plot(t0,ia),grid
ylabel('ia')
xlabel('t')
xlim([0 0.1])

subplot(3,2,4)
plot(t0,iaf),grid
ylabel('ia after LPF')
xlabel('t')
xlim([0 0.1])

subplot(3,2,6)
plot(t,iax),grid
ylabel('ia after COSF')
xlabel('t')
xlim([0 0.1])

figure(6);
subplot(3,2,1)
plot(t0,va),grid;hold on;plot(t,abs(Va),'r')
ylabel('va, Va in volts')
xlabel('Time in seconds')

subplot(3,2,2)
plot(t0,ia),grid;hold on;plot(t,abs(Ia),'r')
ylabel('ia, Ia in amps')
xlabel('Time in seconds')

subplot(3,2,3)
plot(t0,vb),grid;hold on;plot(t,abs(Vb),'r')
ylabel('vb, Vb in volts')
xlabel('Time in seconds')

subplot(3,2,4)
plot(t0,ib),grid;hold on;plot(t,abs(Ib),'r')
ylabel('ib, Ib in amps')
xlabel('Time in seconds')

subplot(3,2,5)
plot(t0,vc),grid;hold on;plot(t,abs(Vc),'r')
ylabel('vc, Vc in volts')
xlabel('Time in seconds')

subplot(3,2,6)
plot(t0,ic),grid;hold on;plot(t,abs(Ic),'r')
ylabel('ic, Ic in amps')
xlabel('Time in seconds')

figure(7);
subplot(3,2,1)
plot(t,angle(Va)*kk,'r')
ylabel('angle(Va) in deg')
xlabel('Time in seconds')

subplot(3,2,2)
plot(t,angle(Ia)*kk,'r')
ylabel('angle(Ia) in deg')
xlabel('Time in seconds')

subplot(3,2,3)
plot(t,angle(Vb)*kk,'r')
ylabel('angle(Vb) in deg')
xlabel('Time in seconds')

subplot(3,2,4)
plot(t,angle(Ib)*kk,'r')
ylabel('angle(Ib) in deg')
xlabel('Time in seconds')

subplot(3,2,5)
plot(t,angle(Vc)*kk,'r')
ylabel('angle(Vc) in deg')
xlabel('Time in seconds')

subplot(3,2,6)
plot(t,angle(Ic)*kk,'r')
ylabel('angle(Ic) in deg')
xlabel('Time in seconds')

figure(8);
subplot(3,2,1)
plot(t,angle(Va)*kk - angle(Va)*kk,'r')
ylabel('angle(Va) in deg')
xlabel('Time in seconds')

subplot(3,2,2)
plot(t,angle(Ia)*kk - angle(Va)*kk,'r')
ylabel('angle(Ia) in deg')
xlabel('Time in seconds')

subplot(3,2,3)
plot(t,angle(Vb)*kk - angle(Va)*kk,'r')
ylabel('angle(Vb) in deg')
xlabel('Time in seconds')

subplot(3,2,4)
plot(t,angle(Ib)*kk - angle(Va)*kk,'r')
ylabel('angle(Ib) in deg')
xlabel('Time in seconds')

subplot(3,2,5)
plot(t,angle(Vc)*kk - angle(Va)*kk,'r')
ylabel('angle(Vc) in deg')
xlabel('Time in seconds')

subplot(3,2,6)
plot(t,angle(Ic)*kk - angle(Va)*kk,'r')
ylabel('angle(Ic) in deg')
xlabel('Time in seconds')
