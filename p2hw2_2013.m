%%% Homework 3 - Protection II - Spring 2012
% Simple script for solution after ATP simulation
% Version for OCTAVE V3.2 and MATLAB 04/1/2012
% by Luis G. Pérez
%
dos('del resulthw3.txt')
diary('resulthw3.txt')
clc
clear all
close all
CTR=800/5;
VTR=525000/115;
%% 1. Read the csv file with the results of the simulation
Ain=csvread('prot2_hw2_2108.csv'); % Ain is the input table
tor=Ain(:,1);  % first column contains time in seconds ("tor" means original time vector)
va=Ain(:,2)/VTR; % Second column contains instantanoues va voltage...
vb=Ain(:,3)/VTR; % Etc...
vc=Ain(:,4)/VTR;
ia=Ain(:,5)/CTR;
ib=Ain(:,6)/CTR;
ic=Ain(:,7)/CTR;
Dtor=tor(2)-tor(1); % Original sampling rate (in seconds)
fsor=1/Dtor;  % Original sampling frequency
%
%% 2. Pass the signals through a simple RC LPF

fc=600;
vaf=rcfilter2(fc,va,tor); % rcfilter is a home-made function to replace lsim
vbf=rcfilter2(fc,vb,tor);
vcf=rcfilter2(fc,vc,tor);
iaf=rcfilter2(fc,ia,tor);
ibf=rcfilter2(fc,ib,tor);
icf=rcfilter2(fc,ic,tor);
%
%% 3. Re-sampling of the waves at 24 samples per cycle
Nsc=16; % Number of samples per cycle with the relay processing algorithm
fs=Nsc*60;  % New sampling frequency in Hz
Dt=1/fs; % New sampling rate in seconds
t=[0:Dt:tor(length(tor))]; % New time vector
vafr=interp1(tor,vaf,t); % Resampling of vaf to the new sampling rate
vbfr=interp1(tor,vbf,t); % Etc
vcfr=interp1(tor,vcf,t);
iafr=interp1(tor,iaf,t);
ibfr=interp1(tor,ibf,t);
icfr=interp1(tor,icf,t);
%
%% 4. Cosine and sine filter: calculation of the filter coeficients and filtering
k=[1:Nsc];
bc=(2/Nsc)*cos(2*pi*k/Nsc);  % Cosine filter coefficients (Nsc samples per cycle)
bs=(2/Nsc)*sin(2*pi*k/Nsc);  % Sine filter coefficients (Nsc samples per cycle)
vax=filter(bc,1,vafr);  % Cosine filtering
vbx=filter(bc,1,vbfr);
vcx=filter(bc,1,vcfr);
vay=filter(bs,1,vafr);  % Sine filtering
vby=filter(bs,1,vbfr);
vcy=filter(bs,1,vcfr);
iax=filter(bc,1,iafr);  % Cosine filtering
ibx=filter(bc,1,ibfr);
icx=filter(bc,1,icfr);
iay=filter(bs,1,iafr);  % Sine filtering
iby=filter(bs,1,ibfr);
icy=filter(bs,1,icfr);
%
%% 5. Static phasor calculation (notice division by sqrt(2))
j=sqrt(-1);
Va=(1/sqrt(2))*(vax+j*vay).*exp(-j*2*pi*60*t);  % Phasor for va
Vb=(1/sqrt(2))*(vbx+j*vby).*exp(-j*2*pi*60*t);
Vc=(1/sqrt(2))*(vcx+j*vcy).*exp(-j*2*pi*60*t);
Ia=(1/sqrt(2))*(iax+j*iay).*exp(-j*2*pi*60*t);  % Phasor for ia
Ib=(1/sqrt(2))*(ibx+j*iby).*exp(-j*2*pi*60*t);
Ic=(1/sqrt(2))*(icx+j*icy).*exp(-j*2*pi*60*t);
%% 6. Figures

figure(1)
subplot(211)
plot(60*tor,va,60*tor,vb,60*tor,vc)
grid on
title('Original instantaneous voltages')
ylabel('Voltage in volts')
xlabel('Time in cycles')
xlim([0 6])
%
subplot(212)
plot(60*tor,ia,60*tor,ib,60*tor,ic)
grid on
title('Original instantaneous currents')
ylabel('Current in amps')
xlabel('Time in cycles')
xlim([0 6])
set(gca,'YMinorGrid','on')
%
figure(2);plot(tor,va),grid;hold on;plot(t,abs(Va),'r')
title('Original instantaneous phase-a voltage and magnitude of the phasor')
ylabel('Voltage in volts')
xlabel('Time in seconds')
%
figure(3);plot(tor,ia),grid;hold on;plot(t,abs(Ia),'r')
title('Original instantaneous phase-a current and magnitude of the phasor')
ylabel('Current in amps')
xlabel('Time in seconds')
%
kk=180/pi;
figure(4);plot(t,angle(Va)*kk,t,angle(Ia)*kk)
title('Angles in degrees: Blue-->Voltage, Green-->Current')
ylabel('Phase-a voltage and current angles in degrees')
xlabel('Time in seconds')
%
figure(5)
subplot(3,2,1)
plot(tor,va),grid
ylabel('va')
xlabel('t')
xlim([0 0.1])

%
subplot(3,2,3)
plot(tor,vaf),grid
ylabel('va after LPF')
xlabel('t')
xlim([0 0.1])
%
subplot(3,2,5)
plot(t,vax),grid
ylabel('va after COSF')
xlabel('t')
xlim([0 0.1])
%
subplot(3,2,2)
plot(tor,ia),grid
ylabel('ia')
xlabel('t')
xlim([0 0.1])
%
subplot(3,2,4)
plot(tor,iaf),grid
ylabel('ia after LPF')
xlabel('t')
xlim([0 0.1])
%
subplot(3,2,6)
plot(t,iax),grid
ylabel('ia after COSF')
xlabel('t')
xlim([0 0.1])
%
figure(6);
subplot(3,2,1)
plot(tor,va),grid;hold on;plot(t,abs(Va),'r')
ylabel('va, Va in volts')
xlabel('Time in seconds')
%
subplot(3,2,2)
plot(tor,ia),grid;hold on;plot(t,abs(Ia),'r')
ylabel('ia, Ia in amps')
xlabel('Time in seconds')
%
subplot(3,2,3)
plot(tor,vb),grid;hold on;plot(t,abs(Vb),'r')
ylabel('vb, Vb in volts')
xlabel('Time in seconds')
%
subplot(3,2,4)
plot(tor,ib),grid;hold on;plot(t,abs(Ib),'r')
ylabel('ib, Ib in amps')
xlabel('Time in seconds')
%
subplot(3,2,5)
plot(tor,vc),grid;hold on;plot(t,abs(Vc),'r')
ylabel('vc, Vc in volts')
xlabel('Time in seconds')
%
subplot(3,2,6)
plot(tor,ic),grid;hold on;plot(t,abs(Ic),'r')
ylabel('ic, Ic in amps')
xlabel('Time in seconds')
%
%
figure(7);
subplot(3,2,1)
plot(t,angle(Va)*kk,'r')
ylabel('angle(Va) in deg')
xlabel('Time in seconds')
%
subplot(3,2,2)
plot(t,angle(Ia)*kk,'r')
ylabel('angle(Ia) in deg')
xlabel('Time in seconds')
%
subplot(3,2,3)
plot(t,angle(Vb)*kk,'r')
ylabel('angle(Vb) in deg')
xlabel('Time in seconds')
%
subplot(3,2,4)
plot(t,angle(Ib)*kk,'r')
ylabel('angle(Ib) in deg')
xlabel('Time in seconds')
%
subplot(3,2,5)
plot(t,angle(Vc)*kk,'r')
ylabel('angle(Vc) in deg')
xlabel('Time in seconds')
%
subplot(3,2,6)
plot(t,angle(Ic)*kk,'r')
ylabel('angle(Ic) in deg')
xlabel('Time in seconds')
%


% Display results in working area
%
kk=180/pi;
disp('Homework 3: Results')
disp('Voltage and current phasors as measured by the relay')
disp(['Va = ' num2str(abs(Va(length(Va)-1))) ' V /angle ' num2str(angle(Va(length(Va)-1))*kk) ' deg'])
disp(['Vb = ' num2str(abs(Vb(length(Vb)-1))) ' V /angle ' num2str(angle(Vb(length(Vb)-1))*kk) ' deg'])
disp(['Vc = ' num2str(abs(Vc(length(Vc)-1))) ' V /angle ' num2str(angle(Vc(length(Vc)-1))*kk) ' deg'])
disp(['Ia = ' num2str(abs(Ia(length(Va)-1))) ' A /angle ' num2str(angle(Ia(length(Va)-1))*kk) ' deg'])
disp(['Ib = ' num2str(abs(Ib(length(Va)-1))) ' A /angle ' num2str(angle(Ib(length(Va)-1))*kk) ' deg'])
disp(['Ic = ' num2str(abs(Ic(length(Va)-1))) ' A /angle ' num2str(angle(Ic(length(Va)-1))*kk) ' deg'])
disp('These results must be compared with the results obtained with')
disp('a traditional fault calculation...')

% Approximate fault results (neglecting the line capacitance)

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
% Use home-made zpolar function to determine the results in polar form
disp(' ')
disp('The results using traditional fault calculations')
disp('and neglecting capacitances are as follows...')
disp(' ')
disp(['Var = ' num2str(abs(Var)) ' V /angle ' num2str(angle(Var)*kk) ' deg'])
disp(['Vbr = ' num2str(abs(Vbr)) ' V /angle ' num2str(angle(Vbr)*kk) ' deg'])
disp(['Vcr = ' num2str(abs(Vcr)) ' V /angle ' num2str(angle(Vcr)*kk) ' deg'])
disp(['Iar = ' num2str(abs(Iar)) ' A /angle ' num2str(angle(Iar)*kk) ' deg'])
disp(['Ibr = ' num2str(abs(Ibr)) ' A /angle ' num2str(angle(Ibr)*kk) ' deg'])
disp(['Icr = ' num2str(abs(Icr)) ' A /angle ' num2str(angle(Icr)*kk) ' deg'])
diary off
