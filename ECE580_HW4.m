%Emanuel Amaro-Zurita
%ECE 580
%Homework 4-5

clc
clear all

CTR=800/5;
VTR=(500000/sqrt(3))/67;
ZTR = VTR/CTR;
pol=180/pi;

zL0 = 0.1+2.6*j;
zL1 = 0.073+j*0.8;
ko=(zL0-zL1)/(3*zL1);
zL1mag = abs(zL1);
zL1ang = pol*angle(zL1);
Zlt = 100*zL1;
Zr = (0.8*zL1*100)/ZTR;
Zrmag = abs(Zr);
Zrang = pol*angle(Zr);
r = 0.5*Zrmag;

%Mho Relay for Zr
theta = 0:pi/50:2*pi;
Zcx = r*cos(theta)+(0.5*real(Zr));
Zcy = r*sin(theta)+(0.5*imag(Zr));
%

%% Problem 1

choice1 = menu('Fault Location','F1','F2');
choice2 = menu('Fault Type','a-g','b-c','b-c-g');
choice3 = menu('Fault Resistance','Do Not Include','Include');
%

if choice1==1 && choice2==1 && choice3 == 1
    fprintf('Case 1\n')
    in = csvread('f1_ag_r0.csv');
end

if choice1==1 && choice2==1 && choice3 == 2
    fprintf('Case 2\n')
    in = csvread('f1_ag_r2.csv');
end

if choice1==1 && choice2==2 && choice3 == 1
    fprintf('Case 3\n')
    in = csvread('f1_bc_r0.csv');
end

if choice1==1 && choice2==2 && choice3 == 2
    fprintf('Case 4\n')
    in = csvread('f1_bc_r2.csv');
end

if choice1==2 && choice2==1 && choice3 == 1
    fprintf('Case 5\n')
    in = csvread('f2_ag_r0.csv');
end

if choice1==2 && choice2==1 && choice3 == 2
    fprintf('Case 6\n')
    in = csvread('f2_ag_r2.csv');
end

if choice1==2 && choice2==2 && choice3 == 1
    fprintf('Case 7\n')
    in = csvread('f2_bc_r0.csv');
end

if choice1==2 && choice2==2 && choice3 == 2
    fprintf('Case 8\n')
    in = csvread('f2_bc_r2.csv');
end

if choice1==1 && choice2==3 && choice3 == 1
    fprintf('Case 9\n')
    in = csvread('f1_bcg_r0.csv');
end

tor = in(:,1);  % original time
va = in(:,2)/VTR; % Instantanoues va Voltage
vb = in(:,3)/VTR; % Instantanoues vb Voltage
vc = in(:,4)/VTR; % Instantanoues vc Voltage
ia = in(:,5)/CTR; % Instantanoues ia Current
ib = in(:,6)/CTR; % Instantanoues ib Current
ic = in(:,7)/CTR; % Instantanoues ic Current
Dtor = tor(2)-tor(1); % Original sampling rate
fsor = 1/Dtor;  % Original sampling frequency55

% Part B: Implement an analog low pass filter with an fc of 600 Hz
fc = 600; % Cutoff Frequency
RC = 1/(2*pi*fc);
Dt = 1E-5;

%Filtered instantanoues va Voltage
vaf(1:length(tor),1) = 0;
for k = 2:length(tor)
vaf(k) = (Dt*va(k)+RC*vaf(k-1))/(Dt+RC);
end

%Filtered instantanoues vb Voltage
vbf(1:length(tor),1) = 0;
for k = 2:length(tor)
vbf(k) = (Dt*vb(k)+RC*vbf(k-1))/(Dt+RC);
end

%Filtered instantanoues vc Voltage
vcf(1:length(tor),1) = 0;
for k=2:length(tor)
vcf(k)=(Dt*vc(k)+RC*vcf(k-1))/(Dt+RC);
end

%Filtered instantanoues ia Current
iaf(1:length(tor),1) = 0;
for k = 2:length(tor)
iaf(k) = (Dt*ia(k)+RC*iaf(k-1))/(Dt+RC);
end

%Filtered instantanoues ib Current
ibf(1:length(tor),1) = 0;
for k = 2:length(tor)
ibf(k) = (Dt*ib(k)+RC*ibf(k-1))/(Dt+RC);
end

%Filtered instantanoues ic Current
icf(1:length(tor),1) = 0;
for k = 2:length(tor)
icf(k) = (Dt*ic(k)+RC*icf(k-1))/(Dt+RC);
end
% Part C: Sample the signal at 16 samples per cycle
Ns = 16;
fs = Ns*60;  % New sampling frequency in Hz
Dts = 1/fs; % New sampling rate in seconds
t = [0:Dts:tor(length(tor))]; % New time vector
vafs = interp1(tor,vaf,t); % Resampling of vaf 
vbfs = interp1(tor,vbf,t); % Resampling of vbf 
vcfs = interp1(tor,vcf,t); % Resampling of vcf 
iafs = interp1(tor,iaf,t); % Resampling of iaf 
ibfs = interp1(tor,ibf,t); % Resampling of ibf 
icfs = interp1(tor,icf,t); % Resampling of icf 
%
% Part D: Implement digital filters in order to obtain the fundamental of the six quantities
m = [1:Ns];
bc = (2/Ns)*cos(2*pi*m/Ns);  % Cosine filter coefficients (16 samples per cycle)
bs = (2/Ns)*sin(2*pi*m/Ns);  % Sine filter coefficients (16 samples per cycle)
vac = filter(bc,1,vafs);  % Cosine filtering
vbc = filter(bc,1,vbfs);
vcc = filter(bc,1,vcfs);
vas = filter(bs,1,vafs);  % Sine filtering
vbs = filter(bs,1,vbfs);
vcs = filter(bs,1,vcfs);
iac = filter(bc,1,iafs);  % Cosine filtering
ibc = filter(bc,1,ibfs);
icc = filter(bc,1,icfs);
ias = filter(bs,1,iafs);  % Sine filtering
ibs = filter(bs,1,ibfs);
ics = filter(bs,1,icfs);

% Part E: Complete your program in a way that it determines the magnitudes and phases of the phasors associated with the six quantities
j=sqrt(-1);
Va=(1/sqrt(2))*(vac+j*vas).*exp(-j*2*pi*60*t);  % Phasor for va
Vb=(1/sqrt(2))*(vbc+j*vbs).*exp(-j*2*pi*60*t);  % Phasor for vb
Vc=(1/sqrt(2))*(vcc+j*vcs).*exp(-j*2*pi*60*t);  % Phasor for vc
Ia=(1/sqrt(2))*(iac+j*ias).*exp(-j*2*pi*60*t);  % Phasor for ia
Ib=(1/sqrt(2))*(ibc+j*ibs).*exp(-j*2*pi*60*t);  % Phasor for ib
Ic=(1/sqrt(2))*(icc+j*ics).*exp(-j*2*pi*60*t);  % Phasor for ic



%Cross Polarizaton
Ires1 = (Ia+Ib+Ic);
Vab1 = Va-Vb;
Vbc1 = Vb-Vc;
Vca1 = Vc-Va;
Iab1 = Ia-Ib;
Ibc1 = Ib-Ic;
Ica1 = Ic-Ia;

Iag1 = (Ia+(ko.*Ires1));
Ibg1 = (Ib+(ko.*Ires1));
Icg1 = (Ic+(ko.*Ires1));

Zab1 = (Vab1)./Iab1;
Zbc1 = (Vbc1)./Ibc1;
Zca1 = (Vca1)./Ica1;
Zag1 = (Va)./(Iag1);
Zbg1 = (Vb)./(Ibg1);
Zcg1 = (Vc)./(Icg1);

S1ab1 = ((Zr.*Iab1)-Vab1);
S1bc1 = ((Zr.*Ibc1)-Vbc1);
S1ca1 = ((Zr.*Ica1)-Vca1);
S1ag1 = ((Zr.*Iag1)-Va);
S1bg1 = ((Zr.*Ibg1)-Vb);
S1cg1 = ((Zr.*Icg1)-Vc);

S2ab1 = (-j.*Vc); 
S2bc1 = (-j.*Va);
S2ca1 = (-j.*Vb);
S2ag1 = (j.*Vbc1);
S2bg1 = (j.*Vca1);
S2cg1 = (j.*Vab1);

AB11 = (S1ab1./(S2ab1));

BC11 = (S1bc1./(S2bc1));

CA11 = (S1ca1./(S2ca1));

AG11 = (S1ag1./(S2ag1));

BG11 = (S1bg1./(S2bg1));

CG11 = (S1cg1./(S2cg1));

for n = 1:193

if -90 < (angle(AB11(n))*pol)&&(angle(AB11(n))*pol) < 90

MHOABXP(n) = 1;

else MHOABXP(n) = 0;

end

if -90 < (angle(BC11(n))*pol)&&(angle(BC11(n))*pol) < 90

MHOBCXP(n) = 1;

else MHOBCXP(n) = 0;

end

if -90 < (angle(CA11(n))*pol)&&(angle(CA11(n))*pol) < 90

MHOCAXP(n) = 1;

else MHOCAXP(n) = 0;

end

if -90 < (angle(AG11(n))*pol)&&(angle(AG11(n))*pol) < 90

MHOAGXP(n) = 1;

else MHOAGXP(n) = 0;

end

if -90 < (angle(BG11(n))*pol)&&(angle(BG11(n))*pol) < 90

MHOBGXP(n) = 1;

else MHOBGXP(n) = 0;

end

if -90 < (angle(CG11(n))*pol)&&(angle(CG11(n))*pol) < 90

MHOCGXP(n) = 1;

else MHOCGXP(n) = 0;

end

end



if choice1==1 && choice2==1 && choice3 == 1
    V = Va;
    I = Iag1;
    S = S2ag1;
end

if choice1==1 && choice2==1 && choice3 == 2
    V = Va;
    I = Iag1;
    S = S2ag1;
end

if choice1==1 && choice2==2 && choice3 == 1
    V = Vbc1;
    I = Ibc1;
    S = S2bc1;
end

if choice1==1 && choice2==2 && choice3 == 2
    V = Vbc1;
    I = Ibc1;
    S = S2bc1;
end

if choice1==2 && choice2==1 && choice3 == 1
    V = Va;
    I = Iag1;
    S = S2ag1;
end

if choice1==2 && choice2==1 && choice3 == 2
    V = Va;
    I = Iag1;
    S = S2ag1;
end

if choice1==2 && choice2==2 && choice3 == 1
    V = Vbc1;
    I = Ibc1;
    S = S2bc1;
end

if choice1==2 && choice2==2 && choice3 == 2
    V = Vbc1;
    I = Ibc1;
    S = S2bc1;
end

if choice1==1 && choice2==3 && choice3 == 1
    V = Vbc1;
    I = Ibc1;
    S = S2bc1;
end
VpXX = V-(S); %V-S
Zp = VpXX/(I); %V-S/I
Zc = .5*((Zr)+Zp);
rp = .5*abs((Zr)-Zp);
th = 0:pi/50:2*pi;

bp = linspace(0,real(Zlt));

bs = bp/ZTR;

b=linspace(0,1) * (Zlt/ZTR);

Rs1 = rp * cos(theta)+ (real(Zc));

Xs1 = rp * sin(theta)+ (imag(Zc));

figure

hold on
plot(bs,imag(b))
plot(Zcx,Zcy);
plot(Rs1,Xs1);
plot(real(Zag1),imag(Zag1),'k');
plot(real(Zbg1),imag(Zbg1),'g');
plot(real(Zcg1),imag(Zcg1),'r');
plot(real(Zab1),imag(Zab1),'b');
plot(real(Zbc1),imag(Zbc1),'m');
plot(real(Zca1),imag(Zca1),'c');
xlim([-10 10])
ylim([-10 10])

title('MHO Characteristic of the Relays');
xlabel('Real (secondary ohms)');
ylabel('Reactive (secondary ohms)');
legend('Line Impedance','Zr (Self-Polarized)','Zp (Cross-Polarized)','Zag','Zbg','Zcg','Zab','Zbc','Zca');
grid on


figure

subplot(3,1,1)
plot(60*tor,va,'r',60*tor,vb,'k',60*tor,vc,'b')
title('Original Instantaneous Voltages')
xlabel('Time (cyc)')
ylabel('Sec. Voltage (V)')
xlim([0 6])

subplot(3,1,2)
plot(60*tor,ia,'r',60*tor,ib,'k',60*tor,ic,'b')
title('Original Instantaneous Currents')
xlabel('Time (cyc)')
ylabel('Sec. Current (A)')
xlim([0 6])

subplot(3,1,3)
plot(t,MHOABXP,t,MHOBCXP,t,MHOCAXP,t,MHOAGXP,t,MHOBGXP,t,MHOCGXP)

title('Relay Operation');
xlabel('Time (s)');
ylabel('Relay Trip');
ylim([-2 2]);
legend('AB','BC','CA','AG','BG','CG');

E = 5E5/sqrt(3);
Zs1 = 1+10*j-(j*20);
Zs0 = 2+30*j-(j*10);
Zs2 = Zs1;
zL1 = 0.073+0.8*j;
zL0 = 0.1+2.6*j;
zL2 = zL1;
I1 = E/(Zs1+zL1*20+Zs2+zL2*20+Zs0+zL0*20);
I0 = I1;
I2 = I1;

V1 = E-Zs1*I1;
V0 = -Zs0*I0;
V2 = -Zs2*I2;
a = 1*exp(j*120/pol);
Iar = (I0+I1+I2)/CTR;
Ibr = (I0+(a^2)*I1+a*I2)/CTR;
Icr = (I0+a*I1+(a^2)*I2)/CTR;
Var = (V0+V1+V2)/VTR;
Vbr = (V0+(a^2)*V1+a*V2)/VTR;
Vcr = (V0+a*V1+(a^2)*V2)/VTR;

%A Matrix
magVa = abs(Var);
magVb = abs(Vbr);
magVc = abs(Vcr);
phVa = pol*angle(Var);
phVb = pol*angle(Vbr);
phVc = pol*angle(Vcr);
magIa = abs(Iar);
magIb = abs(Ibr);
magIc = abs(Icr);
phIa = pol*angle(Iar);
phIb = pol*angle(Ibr);
phIc = pol*angle(Icr);

pol=180/pi;
fprintf('Homework 3: Results\n')
fprintf('Results using phasors measured by the relay\n\n')
fprintf('Subject   Magnitude    Phase\n')
disp(['Va        ' num2str(abs(Va(length(Va)-1))) ' V    ' num2str(angle(Va(length(Va)-1))*pol) ' deg'])
disp(['Vb        ' num2str(abs(Vb(length(Vb)-1))) ' V    ' num2str(angle(Vb(length(Vb)-1))*pol) ' deg'])
disp(['Vc        ' num2str(abs(Vc(length(Vc)-1))) ' V    ' num2str(angle(Vc(length(Vc)-1))*pol) ' deg'])
disp(['Ia        ' num2str(abs(Ia(length(Ia)-1))) ' A    ' num2str(angle(Ia(length(Ia)-1))*pol) ' deg'])
disp(['Ib        ' num2str(abs(Ib(length(Ib)-1))) ' A    ' num2str(angle(Ib(length(Ib)-1))*pol) ' deg'])
disp(['Ic        ' num2str(abs(Ic(length(Ic)-1))) ' A    ' num2str(angle(Ic(length(Ic)-1))*pol) ' deg'])
fprintf('\n')


fprintf('Results using fault calculation methods\n\n')
fprintf('Subject   Magnitude    Phase\n')
fprintf('Va        %4.4f V     %4.4f deg\n',magVa, phVa)
fprintf('Vb        %4.4f V     %4.4f deg\n',magVb, phVb)
fprintf('Vc        %4.4f V     %4.4f deg\n',magVc, phVc)
fprintf('Ia        %4.4f A     %4.4f deg\n',magIa, phIa)
fprintf('Ib        %4.4f A     %4.4f deg\n',magIb, phIb)
fprintf('Ic        %4.4f A     %4.4f deg\n',magIc, phIc)


% 1D: Finding the fault location of F1
figure
m = imag(V./I)./imag(Zlt);
plot(t,m)
title('Plot of m over time');
xlabel('m (per unit percentage)');
ylabel('t (s)');


%% Problem 2
VTR2 = (115000/sqrt(3))/67;
CTR2 = 600/5;

%From 1 to 2
R1_12 = 2.353;
X1_12 = 10.716;

R0_12 = 3.00;
X0_12 = 30.00;

L_12 = 23.3;

zL0_12 = R0_12+X0_12;
zL1_12 = R1_12+X1_12;

ko_12=(zL0_12-zL1_12)/(3*zL1_12);

%From 1 to 3
R1_13 = 0.135;
X1_13 = 1.438;

R0_13 = 0.50;
X0_13 = 5.00;

L_13 = 4.48;

zL0_13 = R0_13+X0_13;
zL1_13 = R1_13+X1_13;

ko_13=(zL0_13-zL1_13)/(3*zL1_13);

%From 3 to 4
R1_34 = 0.098;
X1_34 = 1.039;

R0_34 = 0.20;
X0_34 = 3.50;

L_34 = 3.24;

zL0_34 = R0_34+X0_34;
zL1_34 = R1_34+X1_34;

ko_34=(zL0_34-zL1_34)/(3*zL1_34);

%From 4 to 5
R1_45 = 0.481;
X1_45 = 3.747;

R0_45 = 1.00;
X0_45 = 10.00;

L_45 = 8.00;

zL0_45 = R0_45+X0_45;
zL1_45 = R1_45+X1_45;

ko_45=(zL0_45-zL1_45)/(3*zL1_45);

%From 2 to 5
R1_25 = 1.493;
X1_25 = 6.807;

R0_25 = 2.00;
X0_25 = 17.50;

L_25 = 14.8;

zL0_25 = R0_25+X0_25;
zL1_25 = R1_25+X1_25;

ko_25=(zL0_25-zL1_25)/(3*zL1_25);

%From 1 to 4
R1_14 = 0.238;
X1_14 = 2.477;

R0_14 = 0.8;
X0_14 = 8.00;

L_14 = 7.72;

zL0_14 = R0_14+X0_14;
zL1_14 = R1_14+X1_14;

ko_14=(zL0_14-zL1_14)/(3*zL1_14);

%Zone Calculation no L because units are in Ohms
Zone1_12 = 0.8*(zL1_12);
Zone2_12 = zL1_12*+.5*zL1_25;
fprintf('Problem 2\n')
fprintf('Zone 1 for line 1-2 is %4.4f\n',Zone1_12);
fprintf('Zone 2 for line 1-2 is %4.4f\n',Zone2_12);

Zone1_13 = 0.8*(zL1_13);
Zone2_13 = zL1_13*+.5*zL1_34;
fprintf('\n')
fprintf('Zone 1 for line 1-3 is %4.4f\n',Zone1_13);
fprintf('Zone 2 for line 1-3 is %4.4f\n',Zone2_13);

Zone1_14 = 0.8*(zL1_14);
Zone2_14 = zL1_14+.5*zL1_34;
fprintf('\n')
fprintf('Zone 1 for line 1-4 is %4.4f\n',Zone1_14);
fprintf('Zone 2 for line 1-4 is %4.4f\n',Zone2_14);

Zone1_34 = 0.8*(zL1_34);
Zone2_34 = zL1_34+.5*zL1_45;
fprintf('\n')
fprintf('Zone 1 for line 1-4 is %4.4f\n',Zone1_34);
fprintf('Zone 2 for line 1-4 is %4.4f\n',Zone2_34);

Zone1_45 = 0.8*(zL1_45);
Zone2_45 = zL1_45+.5*zL1_25;
fprintf('\n')
fprintf('Zone 1 for line 1-4 is %4.4f\n',Zone1_45);
fprintf('Zone 2 for line 1-4 is %4.4f\n',Zone2_45);

Zone1_25 = 0.8*(zL1_25);
Zone2_25 = zL1_25+.5*zL1_45;
fprintf('\n')
fprintf('Zone 1 for line 1-4 is %4.4f\n',Zone1_25);
fprintf('Zone 2 for line 1-4 is %4.4f\n',Zone2_25);


fprintf('\n')
fprintf('k0 for Relays between bus 1 and 2 are %4.4f\n',ko_12);
fprintf('k0 for Relays between bus 1 and 3 are %4.4f\n',ko_13);
fprintf('k0 for Relays between bus 3 and 4 are %4.4f\n',ko_34);
fprintf('k0 for Relays between bus 4 and 5 are %4.4f\n',ko_45);
fprintf('k0 for Relays between bus 2 and 5 are %4.4f\n',ko_25);
fprintf('k0 for Relays between bus 1 and 4 are %4.4f\n',ko_14);



%