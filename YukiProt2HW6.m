close all
j=sqrt(-1);
kk = 180/pi;
theta = linspace(0,2*pi);
%Transcribe provided curves for thermal limits and UEL from plot into arrays
Ptherm = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
Qtherm = [-.42 -.43 -.44 -.43 -.41 -.39 -.36 -.3 -.22 -.14 0 .42 .6 .67 .74 .77 .8 .83 .85 .87 .88];
PUEL = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
QUEL = [-.42 -.415 -.4 -.39 -.37 -.33 -.28 -.21 -.16 -.08 0];

%Calculate equation for SSSL
Xd = 1.867;
Xdp = 0.2;
Xt = 0.32;
Xe = 0.1;
Xs = Xt + Xe;
S = 29.6;
VL = 13.8e3;
radss = 1/2*(1/Xs+1/Xd);
Qss = 1/2*(1/Xs-1/Xd);

fprintf('\nRadius of SSSL curve: %f p.u., %f MVA\n', radss, radss*S)
fprintf('Offset of SSSL curve: %f p.u., %f Mvar\n', Qss, Qss*S)

%Plot curves in p.u.
plot(Ptherm,Qtherm)
hold on
plot(PUEL,QUEL)
plot(radss*cos(theta),radss*sin(theta)+Qss)
grid on
xlabel('P(p.u.)')
ylabel('Q(p.u.)')
title('Limit curves and SSSL in per unit')

%Convert to impedance plane
CTR = 1600/5;
VTR = 8000/67;
ZTR = VTR/CTR;
Zb = VL^2/(S*10^6);

radssZ = 1/2*(Xs+Xd)*Zb/ZTR;
QssZ = -1/2*(Xd-Xs)*Zb/ZTR;
Ztherm = zeros(1,length(Ptherm));
fprintf('\nRadius of SSSL curve in impedance plane: %f secondary ohms\n', radssZ)
fprintf('Offset of SSSL curve in impedance plane: %f secondary ohms\n', QssZ)

for k = 1:length(Ztherm)
    Ztherm(k) = 1/(conj(Ptherm(k)+j*Qtherm(k)))*Zb/ZTR;
end
ZUEL = zeros(1,length(PUEL));
for k = 1:length(ZUEL)
    ZUEL(k) = 1/(conj(PUEL(k)+j*QUEL(k)))*Zb/ZTR;
end

figure
plot(real(Ztherm),imag(Ztherm))
hold on
grid on
plot(real(ZUEL),imag(ZUEL))
plot(radssZ*cos(theta),radssZ*sin(theta)+QssZ)
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Limit curves and SSSL in impedance plane')

%Calculate inner and outer MHO circles
MHO1 = (1.1*Xd - Xdp/2)*Zb/ZTR;
radMHO1 = MHO1/2;
XMHO1 = -Xdp/2*Zb/ZTR;
cMHO1 = XMHO1 - radMHO1;
fprintf('Inner MHO radius: %f, offset: %f (secondary ohms)\n',radMHO1,cMHO1)
MHO2 = (1.1*Xd+Xt+Xe)*Zb/ZTR;
radMHO2 = MHO2/2;
XMHO2 = (Xe+Xt)*Zb/ZTR;
cMHO2 = XMHO2 - radMHO2;
fprintf('outer MHO radius: %f, offset: %f (secondary ohms)\n',radMHO2,cMHO2)

%Directional element with MTA 120:
xdir = linspace(-30,30);
ydir = xdir*tan(-30/kk);

%Final plot of all curves
figure
plot(radMHO1*cos(theta),radMHO1*sin(theta)+cMHO1)
hold on
grid on
plot(radMHO2*cos(theta),radMHO2*sin(theta)+cMHO2)
plot(real(Ztherm),imag(Ztherm))
plot(real(ZUEL),imag(ZUEL))
plot(radssZ*cos(theta),radssZ*sin(theta)+QssZ)
plot(xdir,ydir)
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('All generator and relay curves (impedance plane)')

%Algorithm to find optimal tapping for minimal M (problem 5)
tap = 1:0.1:20;
I1 = 4.184;
I2 = 4.348;

results = zeros(length(tap),length(tap));
for k = 1:length(tap)
    for n = 1:length(tap)
        results(k,n) = abs(I1/tap(k)-I2/tap(n))/(I1/tap(k))*100;
    end
end

minres = min(results(:));
[row,col]=find(results==minres);
fprintf('Taps for minimal M: Tap 1 = %g; Tap 2 = %g\n',tap(row(1)),tap(col(1)))