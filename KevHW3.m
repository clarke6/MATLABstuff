% 537 HW 3 Part 2 (I borrowed portions from HW2)
clear all
close all
clc

% Given from Problem 2 (Same for problem 3)
Vnom = 525000;    % Volts (525 kV)
CTR = 800/5;      % According to HW 1, problem 1. Unitless.
VTR = 525000/115; % According to HW 1, problem 1. Unitless.
tdel = 1*10^-5;   % Given in HW 1, Problem 1. Seconds.
tfinal = 0.1;     % Given in HW 1, Problem 1. Seconds.
tfault = 0.016;   % Given in HW 1, Problem 1. Seconds.
f0 = 60;          % Assumed for US power grid. Hz.
w0 = 60*2*pi;     % Converted from f0. rad/s.

% HW2 Part A. 
% Loading in CSV file and plotting data.
data = csvread('prot2_hw2_2108.csv',0,0);% Pulls in ATP data from csv file.
sizing = size(data);
t(sizing,1)=0; % Initializes a time vector with 10,001 points. 
k = 1;        % Initiatlizes a loop counting variable.

while k < (sizing(1)+1);
    t(k)=data(k,1);
    va(k)=data(k,2)./VTR;
    vb(k)=data(k,3)./VTR;
    vc(k)=data(k,4)./VTR;
    ia(k)=data(k,5)./CTR;
    ib(k)=data(k,6)./CTR;
    ic(k)=data(k,7)./CTR;
    k = k + 1;
end
% Plotting Secondary Voltages and Currents

%{
figure(1)                    % Plotting Part A, all on one graph.
plot(t, va)
hold on
plot(t, vb)
plot(t, vc)
plot(t, ia)
plot(t, ib)
plot(t, ic)
hold off

figure(2)   % Plotting part A, just the voltage.
plot(t, va)
hold on
plot(t, vb)
plot(t, vc)
hold off

figure(3)   % Plotting part A, just the current.
plot(t, ia)
hold on
plot(t, ib)
plot(t, ic)
hold off
%}

% HW2 Part B. 
% Implementing an analog low-pass RC filter.
fc=600; % Hz
wc=2*pi*fc; % rad/s
RC=1/wc;

% I wrote a function called AnaFilt to perform this recursive formula.
vaAF = AnaFilt(t,va,RC); % Secondary voltage after Analog Filtration.
vbAF = AnaFilt(t,vb,RC);
vcAF = AnaFilt(t,vc,RC);
iaAF = AnaFilt(t,ia,RC); % Secondary current after Analog Filtration.
ibAF = AnaFilt(t,ib,RC);
icAF = AnaFilt(t,ic,RC);

i=sqrt(-1);
%{
figure(4)
subplot(3,2,1)
plot(t,vaAF)
xlabel('Time (s)')
ylabel('Va after LPF (V)')
subplot(3,2,2)
plot(t,iaAF)
xlabel('Time (s)')
ylabel('Ia after LPF (A)')
subplot(3,2,3)
plot(t,vbAF)
xlabel('Time (s)')
ylabel('Vb after LPF (V)')
subplot(3,2,4)
plot(t,ibAF)
xlabel('Time (s)')
ylabel('Ib after LPF (A)')
subplot(3,2,5)
plot(t,vcAF)
xlabel('Time (s)')
ylabel('Vc after LPF (V)')
subplot(3,2,6)
plot(t,icAF)
xlabel('Time (s)')
ylabel('Ic after LPF (A)')
%}



% HW2 Part C.
% Sample the signals at 16 samples per cycles using interpolation.
N = 16;      % This is the number of samples per cycle I will measure.
fs = N * 60; % This is the sampling frequency in Hz.
tdelsamp = 1/fs; % Change in time between samples.
tsamp = [0:tdelsamp:t(sizing)]; % Sampling time vector. New change in time
                              % and new length

vaAFsamp = interp1(t,vaAF,tsamp); % Sampling the filtered secondary voltage
vbAFsamp = interp1(t,vbAF,tsamp); % and interpolating between points.
vcAFsamp = interp1(t,vcAF,tsamp);
iaAFsamp = interp1(t,iaAF,tsamp); % Sampling the filtered secondary current
ibAFsamp = interp1(t,ibAF,tsamp); % and interpolating between points.
icAFsamp = interp1(t,icAF,tsamp);


% HW2 Part D. 
% 1-cycle cosine and sine filters.
q = [1:N]; % Defining q as each integer value of samples per second.
bc = (2/N) * cos(2*pi*q/N);    % Calculating Xq from the cosine formula.
bs = (2/N) * sin(2*pi*q/N);    % Calculating jYq from the sine formula.

vaAFDC = filter(bc,1,vaAFsamp);  % Applying the cosine filter for Xq for v.
vaAFDS = filter(bs,1,vaAFsamp);  % Applying the sine filter for jYq for v.
vbAFDC = filter(bc,1,vbAFsamp);  
vbAFDS = filter(bs,1,vbAFsamp);
vcAFDC = filter(bc,1,vcAFsamp);
vcAFDS = filter(bs,1,vcAFsamp);

iaAFDC = filter(bc,1,iaAFsamp);  % Applying the cosine filter for Xq for i.
iaAFDS = filter(bs,1,iaAFsamp);  % Applying the sine filter for jYq for i.
ibAFDC = filter(bc,1,ibAFsamp);
ibAFDS = filter(bs,1,ibAFsamp);
icAFDC = filter(bc,1,icAFsamp);
icAFDS = filter(bs,1,icAFsamp);

%{
figure(5)
subplot(3,2,1)
plot(tsamp,vaAFDC)
xlabel('Time (s)')
ylabel('Va after COSF (V)')

subplot(3,2,2)
plot(tsamp,iaAFDC)
xlabel('Time (s)')
ylabel('Ia after COSF (A)')

subplot(3,2,3)
plot(tsamp,vbAFDC)
xlabel('Time (s)')
ylabel('Vb after COSF (V)')

subplot(3,2,4)
plot(tsamp,ibAFDC)
xlabel('Time (s)')
ylabel('Ib after COSF (A)')

subplot(3,2,5)
plot(tsamp,vcAFDC)
xlabel('Time (s)')
ylabel('Vc after COSF (V)')

subplot(3,2,6)
plot(tsamp,icAFDC)
xlabel('Time (s)')
ylabel('Ic after COSF (A)')

figure(6)
subplot(3,2,1)
plot(tsamp,vaAFDS)
xlabel('Time (s)')
ylabel('Va after SINF (V)')

subplot(3,2,2)
plot(tsamp,iaAFDS)
xlabel('Time (s)')
ylabel('Ia after SINF (A)')

subplot(3,2,3)
plot(tsamp,vbAFDS)
xlabel('Time (s)')
ylabel('Vb after SINF (V)')

subplot(3,2,4)
plot(tsamp,ibAFDS)
xlabel('Time (s)')
ylabel('Ib after SINF (A)')

subplot(3,2,5)
plot(tsamp,vcAFDS)
xlabel('Time (s)')
ylabel('Vc after SINF (V)')

subplot(3,2,6)
plot(tsamp,icAFDS)
xlabel('Time (s)')
ylabel('Ic after COSF (A)')
%}

% HW2 Part E. 
% Calculate the phasors for va, vb, vc, ia, ib, and ic.
% Note I will calculate the raw values here, and then reference them to a.
vaAFDF = (vaAFDC + j*vaAFDS) .* exp(-j*w0*tsamp) / sqrt(2);  % Complex value for voltage,
vbAFDF = (vbAFDC + j*vbAFDS) .* exp(-j*w0*tsamp) / sqrt(2);  % unreferenced.
vcAFDF = (vcAFDC + j*vcAFDS) .* exp(-j*w0*tsamp) / sqrt(2); 
iaAFDF = (iaAFDC + j*iaAFDS) .* exp(-j*w0*tsamp) / sqrt(2);  % Complex value for current,
ibAFDF = (ibAFDC + j*ibAFDS) .* exp(-j*w0*tsamp) / sqrt(2);  % unreferenced.
icAFDF = (icAFDC + j*icAFDS) .* exp(-j*w0*tsamp) / sqrt(2);

% Converting calculations to magnitudes and phasor angles.
vaAFDFmag = abs(vaAFDF); % Caclulates the magnitude of the corrected va.
vaAFDFang = angle(vaAFDF) .* 57.2958; % Calculates va 
                  %angle. Note this angle will be referenced to 0 degrees.
vbAFDFmag = abs(vbAFDF); % Magnitude of corrected vb.
vbAFDFang = atan2(imag(vbAFDF),real(vbAFDF)) .* 57.2958;
vcAFDFmag = abs(vcAFDF); % Magnitude of corrected vb.
vcAFDFang = atan2(imag(vcAFDF),real(vcAFDF)) .* 57.2958;
iaAFDFmag = abs(iaAFDF);
iaAFDFang = atan2(imag(iaAFDF),real(iaAFDF)) .* 57.2958;
ibAFDFmag = abs(ibAFDF);
ibAFDFang = atan2(imag(ibAFDF),real(ibAFDF)) .* 57.2958;
icAFDFmag = abs(icAFDF);
icAFDFang = atan2(imag(icAFDF),real(icAFDF)) .* 57.2958;

% Rereferencing phasors to the va phase angle.

vaAFDFangRR = vaAFDFang - vaAFDFang;
vbAFDFangRR = vbAFDFang - vaAFDFang;
vcAFDFangRR = vcAFDFang - vaAFDFang;
iaAFDFangRR = iaAFDFang - vaAFDFang;
ibAFDFangRR = ibAFDFang - vaAFDFang;
icAFDFangRR = icAFDFang - vaAFDFang;


%{
% Additional Plots
figure(7);plot(t,va),grid;hold on;plot(tsamp,vaAFDFmag,'r')
title('Original instantaneous phase-a voltage and magnitude of the filtered phasor')
ylabel('Voltage in volts')
xlabel('Time in seconds')

figure(8);plot(t,ia),grid;hold on;plot(tsamp,iaAFDFmag,'r')
title('Original instantaneous phase-a current and magnitude of the phasor')
ylabel('Current in amps')
xlabel('Time in seconds')

figure(9) 
plot(tsamp,vaAFDFang) 
hold on 
plot(tsamp, iaAFDFang)
title('Angles in degrees: Blue-->Voltage, Green-->Current')
ylabel('Phase-a voltage and current angles in degrees')
xlabel('Time in seconds')


figure(10)
subplot(3,2,1)
plot(t,va),grid;hold on;plot(tsamp,vaAFDFmag,'r')
ylabel('va, Va in volts')
xlabel('Time in seconds')

subplot(3,2,2)
plot(t,ia),grid;hold on;plot(tsamp,iaAFDFmag,'r')
ylabel('ia, Ia in amps')
xlabel('Time in seconds')

subplot(3,2,3)
plot(t,vb),grid;hold on;plot(tsamp,vbAFDFmag,'r')
ylabel('vb, Vb in volts')
xlabel('Time in seconds')

subplot(3,2,4)
plot(t,ib),grid;hold on;plot(tsamp,ibAFDFmag,'r')
ylabel('ib, Ib in amps')
xlabel('Time in seconds')

subplot(3,2,5)
plot(t,vc),grid;hold on;plot(tsamp,vcAFDFmag,'r')
ylabel('vc, Vc in volts')
xlabel('Time in seconds')

subplot(3,2,6)
plot(t,ic),grid;hold on;plot(tsamp,icAFDFmag,'r')
ylabel('ic, Ic in amps')
xlabel('Time in seconds')


% Referencing va to zero.
figure(11);
subplot(3,2,1)
plot(tsamp,vaAFDFangRR,'r')
ylabel('angle(Va) in deg')
xlabel('Time in seconds')

subplot(3,2,2)
plot(tsamp,iaAFDFangRR,'r')
ylabel('angle(Ia) in deg')
xlabel('Time in seconds')

subplot(3,2,3)
plot(tsamp,vbAFDFangRR,'r')
ylabel('angle(Vb) in deg')
xlabel('Time in seconds')

subplot(3,2,4)
plot(tsamp,ibAFDFangRR,'r')
ylabel('angle(Ib) in deg')
xlabel('Time in seconds')

subplot(3,2,5)
plot(tsamp,vcAFDFangRR,'r')
ylabel('angle(Vc) in deg')
xlabel('Time in seconds')

subplot(3,2,6)
plot(tsamp,icAFDFangRR,'r')
ylabel('angle(Ic) in deg')
xlabel('Time in seconds')

% Without referencing vA to zero
figure(12);
subplot(3,2,1)
plot(tsamp,vaAFDFang,'r')
ylabel('angle(Va) in deg')
xlabel('Time in seconds')

subplot(3,2,2)
plot(tsamp,iaAFDFang,'r')
ylabel('angle(Ia) in deg')
xlabel('Time in seconds')

subplot(3,2,3)
plot(tsamp,vbAFDFang,'r')
ylabel('angle(Vb) in deg')
xlabel('Time in seconds')

subplot(3,2,4)
plot(tsamp,ibAFDFang,'r')
ylabel('angle(Ib) in deg')
xlabel('Time in seconds')

subplot(3,2,5)
plot(tsamp,vcAFDFang,'r')
ylabel('angle(Vc) in deg')
xlabel('Time in seconds')

subplot(3,2,6)
plot(tsamp,icAFDFang,'r')
ylabel('angle(Ic) in deg')
xlabel('Time in seconds')
%}

% Calculating the fault using conventional phasor methods. Again using
% the a phase of voltage as aa reference.
E = 525000 / sqrt(3);
Zs1 = 2+30*j;
Zs2 = Zs1;
Zs0 = Zs1;
zL1 = 0.073+0.8*j;
zL2 = zL1;
zL0 = 0.1+2.3*j;
len = 20; % Miles from the fault.
ZL1 = len * zL1;
ZL2 = len * zL2;
ZL0 = len * zL0;

% Calculating Sequence Components for V and I
I1 = E / (Zs1 + ZL1 + Zs2 + ZL2 + Zs0 + ZL0);
I2 = I1;
I0 = I1;
V1 = E - Zs1 * I1;  % Positive sequence has E included.
V2= -Zs1 * I2;
V0= -Zs0 * I0;

a= 1 * exp(j*120/57.2958);
Iar = (I0 + I1 + I2) / CTR;
Ibr = (I0 + (a^2)*I1 + a*I2) / CTR;
Icr = (I0 + a*I1 + (a^2)*I2) / CTR;
Var = (V0 + V1 + V2) / VTR;
Vbr = (V0 + (a^2)*V1 + a*V2) /VTR;
Vcr = (V0 + a*V1 + (a^2)*V2) /VTR;

% HW3, Problem 1. Caluclating MHO relay inputs.
% Use home-made zpolar function to determine the results in polar form
disp(' ')
disp('Homework 3, Problem 1: Results')
disp('The secondary phase voltages and currents results using traditional ')
disp('fault calculations and neglecting capacitances are as follows...')
disp(' ')
disp(['Var = ' num2str(abs(Var)) ' V /angle ' num2str(angle(Var)*57.2958) ' deg'])
disp(['Vbr = ' num2str(abs(Vbr)) ' V /angle ' num2str(angle(Vbr)*57.2958) ' deg'])
disp(['Vcr = ' num2str(abs(Vcr)) ' V /angle ' num2str(angle(Vcr)*57.2958) ' deg'])
disp(['Iar = ' num2str(abs(Iar)) ' A /angle ' num2str(angle(Iar)*57.2958) ' deg'])
disp(['Ibr = ' num2str(abs(Ibr)) ' A /angle ' num2str(angle(Ibr)*57.2958) ' deg'])
disp(['Icr = ' num2str(abs(Icr)) ' A /angle ' num2str(angle(Icr)*57.2958) ' deg'])


% Phase Distance Elements
Vab = Var - Vbr;
Vbc = Vbr - Vcr;
Vca = Vcr - Var;
Iab = Iar - Ibr;
Ibc = Ibr - Icr;
Ica = Icr - Iar;

% Ground Distance Elements
k0 = (ZL0 - ZL1) / (3*ZL1);
Ires = Iar + Ibr + Icr;
Vag = Var;
Vbg = Vbr;
Vcg = Vcr;
Iag = Iar + k0*Ires;
Ibg = Ibr + k0*Ires;
Icg = Icr + k0*Ires;

% Apparent Impedances of each element.
Zab = Vab/Iab;
Zbc = Vbc/Ibc;
Zca = Vca/Ica;
Zag = Vag/Iag;
Zbg = Vbg/Ibg;
Zcg = Vcg/Icg;


disp(' ')
disp('The secondary apparent impdedances using traditional methods are ')
disp('as follows ')
disp(' ')
disp(['Zab = ' num2str(abs(Zab)) ' Ohms Secondary / angle ' num2str(angle(Zab)*57.2958) ' deg / Rab = ' num2str(real(Zab)) ' ohms secondary / Xab = ' num2str(imag(Zab)) ' ohms secondary'])
disp(['Zbc = ' num2str(abs(Zbc)) ' Ohms Secondary / angle ' num2str(angle(Zbc)*57.2958) ' deg / Rbc = ' num2str(real(Zbc)) ' ohms secondary / Xbc = ' num2str(imag(Zbc)) ' ohms secondary'])
disp(['Zca = ' num2str(abs(Zca)) ' Ohms Secondary / angle ' num2str(angle(Zca)*57.2958) ' deg / Rca = ' num2str(real(Zca)) ' ohms secondary / Xca = ' num2str(imag(Zca)) ' ohms secondary'])
disp(['Zag = ' num2str(abs(Zag)) ' Ohms Secondary / angle ' num2str(angle(Zag)*57.2958) ' deg / Rag = ' num2str(real(Zag)) ' ohms secondary / Xag = ' num2str(imag(Zag)) ' ohms secondary'])
disp(['Zbg = ' num2str(abs(Zbg)) ' Ohms Secondary / angle ' num2str(angle(Zbg)*57.2958) ' deg / Rbg = ' num2str(real(Zbg)) ' ohms secondary / Xbg = ' num2str(imag(Zbg)) ' ohms secondary'])
disp(['Zcg = ' num2str(abs(Zcg)) ' Ohms Secondary / angle ' num2str(angle(Zcg)*57.2958) ' deg / Rbg = ' num2str(real(Zcg)) ' ohms secondary / Xbg = ' num2str(imag(Zcg)) ' ohms secondary'])
disp(' ')

% Since I'm using a self-polarizing scheme, Zp = 0, and then the operation
% criteria becomes -90 <= angle(Zr-Z/Z) <= 90 degrees. For simplicity I
% will simply use < in place of < or = signs.

% Defining Zr, as defined above, zL1 = 0.073 + j*0.8 in unites Ohms/mile
ZTR = VTR/CTR; % For converting from primary to secondary ohms.
L = 100;       % Given, miles
ZL = L * zL1;  % Under normal conditions, only zL1 is present.
ZR = 0.8 * ZL; % ZR triggers at 0.8 of ZL.
ZP = 0;        % Self-polarizing
%ZC = ZR/2
Zr = ZR / ZTR; % Converting from primary to secondary ohms.

% Grpahing
theta=0:0.01:2*pi;       % Defines a range of numbers from 0 to 2pi.
Zrcos=abs(Zr)/2*cos(theta)+real(Zr)/2; % Defines the x values of my circle.
Zrsin=abs(Zr)/2*sin(theta)+imag(Zr)/2; % Defines the y values of my circle.

figure(14)
plot(Zrcos,Zrsin)                % MHO circle with radius Zr.
hold on
plot(real(Zab),imag(Zab),'x')    % Zab apparent impedance
plot(real(Zbc),imag(Zbc),'+')    % Zbc apparent impdedance
plot(real(Zca),imag(Zca),'o')    % Zca
plot(real(Zag),imag(Zag),'*')    % Zag
plot(real(Zbg),imag(Zbg),'s')    % Zbg
plot(real(Zcg),imag(Zcg),'d')    % Zcg


% HW 3, Problem 2. I will be borrowing the k0 and the Zr calculated from
% previous calculations.

disp('')
disp('')
disp('')
disp('Homework 3, Problem 2: Results')
disp('Voltage and current phasors as measured by the relay')
disp(['Va = ' num2str((vaAFDFmag(length(vaAFDFmag)))) ' V /angle ' num2str(vaAFDFang(length(vaAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp(['Vb = ' num2str((vbAFDFmag(length(vbAFDFmag)))) ' V /angle ' num2str(vbAFDFang(length(vbAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp(['Vc = ' num2str((vcAFDFmag(length(vcAFDFmag)))) ' V /angle ' num2str(vcAFDFang(length(vcAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp(['Ia = ' num2str((iaAFDFmag(length(iaAFDFmag)))) ' V /angle ' num2str(iaAFDFang(length(iaAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp(['Ib = ' num2str((ibAFDFmag(length(ibAFDFmag)))) ' V /angle ' num2str(ibAFDFang(length(ibAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp(['Ic = ' num2str((icAFDFmag(length(icAFDFmag)))) ' V /angle ' num2str(icAFDFang(length(icAFDFangRR))-vaAFDFang(length(vaAFDFangRR))) ' deg'])
disp('')


% I will now take the unreferenced values from HW 2, vaAFDF, vbAFDF,
% vcAFDF, iaAFDF, ibAFDF, icAFDF, and again calculate my apparent
% apparent impedances. While for part 1, I solved this problem using
% scalars, I will now have vectors which are changing with time. The
% change in time here is equal to tsamp.

% Digital Phase Distance Elements
VabD = vaAFDF - vbAFDF;
VbcD = vbAFDF - vcAFDF;
VcaD = vcAFDF - vaAFDF;
IabD = iaAFDF - ibAFDF;
IbcD = ibAFDF - icAFDF;
IcaD = icAFDF - iaAFDF;

% Digital Ground Distance Elements (I used calculated k0 for the currents).
IresD = iaAFDF + ibAFDF + icAFDF;
VagD = vaAFDF;
VbgD = vbAFDF;
VcgD = vcAFDF;
IagD = IabD + k0 .* IresD;
IbgD = IbcD + k0 .* IresD;
IcgD = IcaD + k0 .* IresD;

% Digital Apparent Impedances of each element.
ZabD = VabD ./ IabD;
ZbcD = VbcD ./ IbcD;
ZcaD = VcaD ./ IcaD;
ZagD = VagD ./ IagD;
ZbgD = VbgD ./ IbgD;
ZcgD = VcgD ./ IcgD;

% Per HW3, Problem 2 instuctions, I am making binary variables to 
% represent whether a relay operates or not. These will be initially set
% to zero (no-trip). They will be represented as vectors, so that I can
% determine at what time they will trip.
MHOAB = 0 .* ZabD;      % This just ensures I get the right size vectors.
MHOBC = 0 .* ZbcD;
MHOCA = 0 .* ZcaD;
MHOAG = 0 .* ZagD;
MHOBG = 0 .* ZbgD;
MHOCG = 0 .* ZcgD;



k = 1;
kend = length(MHOAB);
%real( (-VabD(k) + IabD(k) * Zr) * conj(VabD(k)))

while (k < kend)
    if ( real( (-VabD(k) + IabD(k) * Zr) * conj(VabD(k)) ) > 0 )
        MHOAB(k) = 1;
    end
    if ( real( (-VbcD(k) + IbcD(k) * Zr) * conj(VbcD(k)) ) > 0 )
        MHOBC(k) = 1;
    end
    if ( real( (-VcaD(k) + IcaD(k) * Zr) * conj(VcaD(k)) ) > 0 )
        MHOCA(k) = 1;
    end
    if ( real( (-VagD(k) + IagD(k) * Zr) * conj(VagD(k)) ) > 0 )
        MHOAG(k) = 1;
    end
    if ( real( (-VbgD(k) + IbgD(k) * Zr) * conj(VbgD(k)) ) > 0 )
        MHOBG(k) = 1;
    end
    if ( real( (-VcgD(k) + IcgD(k) * Zr) * conj(VcgD(k)) ) > 0 )
        MHOCG(k) = 1;
    end
    k = k +1;
end

figure(15)
subplot(3,2,1)
plot(tsamp,MHOAB),grid;hold on;
ylabel('MHOAB: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')

subplot(3,2,2)

plot(tsamp,MHOAG),grid;
ylabel('MHOAG: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')

subplot(3,2,3)
plot(tsamp,MHOBC),grid;hold on;
ylabel('MHOBC: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')


subplot(3,2,4)
plot(tsamp,MHOBG),grid;
ylabel('MHOBG: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')


subplot(3,2,5)
plot(tsamp,MHOCA),grid;hold on;
ylabel('MHOCA: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')


subplot(3,2,6)
plot(tsamp,MHOCG),grid;
ylabel('MHOCG: 0 = No-TRIP, 1 = TRIP')
xlabel('Time in seconds')


figure(16)
plot(Zrcos,Zrsin)                % MHO circle with radius Zr.
hold on
plot(real(ZabD(length(ZabD)-1)),imag(ZabD(length(ZabD)-1)),'x')    % Zab apparent impedance
%plot(real(ZbcD(length(ZbcD)-1)),imag(ZbcD(length(ZbcD)-1)),'+')    % Zbc apparent impdedance
plot(real(ZcaD(length(ZcaD)-1)),imag(ZcaD(length(ZcaD)-1)),'o')    % Zca
plot(real(ZagD(length(ZagD)-1)),imag(ZagD(length(ZagD)-1)),'*')    % Zag
%plot(real(ZbgD(length(ZbgD)-1)),imag(ZbgD(length(ZbgD)-1)),'s')    % Zbg
plot(real(ZcgD(length(ZcgD)-1)),imag(ZcgD(length(ZabD)-1)),'d')    % Zcg

figure(17)
plot(Zrcos,Zrsin)                % MHO circle with radius Zr.
plot(real(ZagD),imag(ZagD))    % Zag

disp(' ')
disp('The secondary apparent impdedances using traditional methods are ')
disp('as follows ')
disp(' ')
disp(['Zab = ' num2str(abs(ZabD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZabD(length(ZabD)-1))*57.2958) ' deg / Rab = ' num2str(real(ZabD(length(ZabD)-1))) ' ohms secondary / Xab = ' num2str(imag(ZabD(length(ZabD)-1))) ' ohms secondary'])
disp(['Zbc = ' num2str(abs(ZbcD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZbcD(length(ZabD)-1))*57.2958) ' deg / Rbc = ' num2str(real(ZbcD(length(ZabD)-1))) ' ohms secondary / Xbc = ' num2str(imag(ZbcD(length(ZabD)-1))) ' ohms secondary'])
disp(['Zca = ' num2str(abs(ZcaD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZcaD(length(ZabD)-1))*57.2958) ' deg / Rca = ' num2str(real(ZcaD(length(ZabD)-1))) ' ohms secondary / Xca = ' num2str(imag(ZcaD(length(ZabD)-1))) ' ohms secondary'])
disp(['Zag = ' num2str(abs(ZagD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZagD(length(ZabD)-1))*57.2958) ' deg / Rag = ' num2str(real(ZagD(length(ZabD)-1))) ' ohms secondary / Xag = ' num2str(imag(ZagD(length(ZabD)-1))) ' ohms secondary'])
disp(['Zbg = ' num2str(abs(ZbgD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZbgD(length(ZabD)-1))*57.2958) ' deg / Rbg = ' num2str(real(ZbgD(length(ZabD)-1))) ' ohms secondary / Xbg = ' num2str(imag(ZbgD(length(ZabD)-1))) ' ohms secondary'])
disp(['Zcg = ' num2str(abs(ZbcD(length(ZabD)-1))) ' Ohms Secondary / angle ' num2str(angle(ZcgD(length(ZabD)-1))*57.2958) ' deg / Rbg = ' num2str(real(ZcgD(length(ZabD)-1))) ' ohms secondary / Xbg = ' num2str(imag(ZcgD(length(ZabD)-1))) ' ohms secondary'])
disp(' ')