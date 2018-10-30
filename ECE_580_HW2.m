%% a data
CTR=800/5;
VTR=(500000/sqrt(3))/67;

Ain=csvread('HW2_Prob1_Fault_Case.csv'); % Ain is the input table
tor=Ain(:,1);  % first column contains time in seconds ("tor" means original time vector)
va=Ain(:,2)/VTR; % Second column contains instantanoues va voltage...
vb=Ain(:,3)/VTR; % Etc...
vc=Ain(:,4)/VTR;
ia=Ain(:,5)/CTR;
ib=Ain(:,6)/CTR;
ic=Ain(:,7)/CTR;
Dtor=tor(2)-tor(1); % Original sampling rate (in seconds)
fsor=1/Dtor;

%% a plot
subplot(211)
plot(tor,va,tor,vb,tor,vc)
grid on
title('3 phase voltages in secondary values')
ylabel('Voltage')
xlabel('Time in seconds')
xlim([0 0.1])

subplot(212)
plot(tor,ia,tor,ib,tor,ic)
grid on
title('3 phase currents in secondary values')
ylabel('Current')
xlabel('Time in seconds')
xlim([0 0.1])

%% b
%va vb vc
fc=600;
RC=1/(2*pi*fc);
Dt=1*(10^-5);  

vaf(1:length(tor),1)=0;
    for k=2:length(tor)
        vaf(k)=(Dt*va(k)+RC*vaf(k-1))/(Dt+RC);
    end
 %plot(tor,vaf)

vbf(1:length(tor),1)=0;
    for k=2:length(tor)
        vbf(k)=(Dt*vb(k)+RC*vbf(k-1))/(Dt+RC);
    end
 %plot(tor,vbf)
 
vcf(1:length(tor),1)=0;
    for k=2:length(tor)
        vcf(k)=(Dt*vc(k)+RC*vcf(k-1))/(Dt+RC);
    end
 %plot(tor,vcf)
  
 %filtered plot V
subplot(211)
plot(tor,vaf,tor,vbf,tor,vcf)
grid on
title('Filtered 3 phase voltages in secondary values with fc=600Hz')
ylabel('Voltage')
xlabel('Time in seconds')
xlim([0 0.1])


 %ia ib ic 
iaf(1:length(tor),1)=0;
    for k=2:length(tor)
        iaf(k)=(Dt*ia(k)+RC*iaf(k-1))/(Dt+RC);
    end
 %plot(tor,iaf)

ibf(1:length(tor),1)=0;
    for k=2:length(tor)
        ibf(k)=(Dt*ib(k)+RC*ibf(k-1))/(Dt+RC);
    end
 %plot(tor,ibf)
 
icf(1:length(tor),1)=0;
    for k=2:length(tor)
        icf(k)=(Dt*ic(k)+RC*icf(k-1))/(Dt+RC);
    end
 %plot(tor,icf)

 %filtered plot I
subplot(212)
plot(tor,iaf,tor,ibf,tor,icf)
grid on
title('Filtered 3 phase current in secondary values with fc=600Hz')
ylabel('Current')
xlabel('Time in seconds')
xlim([0 0.1])

%% c
%3. Re-sampling of the waves at 16 samples per cycle
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

tc=t.*60;
%16 samples per cycle plot V
subplot(211)
plot(tc,vafr,tc,vbfr,tc,vcfr)
grid on
title('Wave of 3 phase voltage with 16 samples per cycle')
ylabel('Voltage')
xlabel('Time in cycles')
%xlim([0 0.1])

%16 samples per cycle plot I
subplot(212)
plot(tc,iafr,tc,ibfr,tc,icfr)
grid on
title('Wave of 3 phase current with 16 samples per cycle')
ylabel('Current')
xlabel('Time in cycles')
%xlim([0 0.1])

%% d
Nsc=16;
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

subplot(221)
plot(tc,vax,tc,vbx,tc,vcx)
grid on
title('Voltage Cosine filter coefficients (16 samples per cycle)')
ylabel('Voltage')
xlabel('Time in cycles')

subplot(222)
plot(tc,vay,tc,vby,tc,vcy)
grid on
title('Voltage Sine filter coefficients (16 samples per cycle)')
ylabel('Voltage')
xlabel('Time in cycles')

subplot(223)
plot(tc,iax,tc,ibx,tc,icx)
grid on
title('Current Cosine filter coefficients (16 samples per cycle)')
ylabel('Current')
xlabel('Time in cycles')

subplot(224)
plot(tc,iay,tc,iby,tc,icy)
grid on
title('Current Sine filter coefficients (16 samples per cycle)')
ylabel('Current')
xlabel('Time in cycles')

%% e
j=sqrt(-1);
Va=(1/sqrt(2))*(vax+j*vay).*exp(-j*2*pi*60*t);  % Phasor for va
Vb=(1/sqrt(2))*(vbx+j*vby).*exp(-j*2*pi*60*t);
Vc=(1/sqrt(2))*(vcx+j*vcy).*exp(-j*2*pi*60*t);
Ia=(1/sqrt(2))*(iax+j*iay).*exp(-j*2*pi*60*t);  % Phasor for ia
Ib=(1/sqrt(2))*(ibx+j*iby).*exp(-j*2*pi*60*t);
Ic=(1/sqrt(2))*(icx+j*icy).*exp(-j*2*pi*60*t);

% Va
% Vb
% Vc
% Ia
% Ib
% Ic
kk=180/pi;
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
disp('Homework 2 Part e: Results')
disp('Voltage and current phasors as measured by the relay')
disp(['Va = ' num2str(abs(Va(length(Va)-1))) ' V /angle ' num2str(angle(Va(length(Va)-1))*kk) ' deg'])
disp(['Vb = ' num2str(abs(Vb(length(Vb)-1))) ' V /angle ' num2str(angle(Vb(length(Vb)-1))*kk) ' deg'])
disp(['Vc = ' num2str(abs(Vc(length(Vc)-1))) ' V /angle ' num2str(angle(Vc(length(Vc)-1))*kk) ' deg'])
disp(['Ia = ' num2str(abs(Ia(length(Va)-1))) ' A /angle ' num2str(angle(Ia(length(Va)-1))*kk) ' deg'])
disp(['Ib = ' num2str(abs(Ib(length(Va)-1))) ' A /angle ' num2str(angle(Ib(length(Va)-1))*kk) ' deg'])
disp(['Ic = ' num2str(abs(Ic(length(Va)-1))) ' A /angle ' num2str(angle(Ic(length(Va)-1))*kk) ' deg'])
disp('These results must be compared with the results obtained with')
disp('a traditional fault calculation...')

%% f
E=500000/sqrt(3);
Zs1=2+30*j;
Zs0=Zs1;
zL1=0.073+0.8*j;
zL0=0.1+2.6*j;
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
