%This program plots the induced torque of an induction machine at two
%specific speeds as a function of the added resistance in a set of three
%ganged potentiometers in the rotor circuit. It also plots the power
%dissipates in the individual potentiometers, so that the necessary
%resistance values and power ratings can be found using the data cursor. It
%then plots four torque-speed curves using evenly spaced potentiometer
%resistance values across the desired range, showing the effect of added
%rotor resistance on the machine's torque-speed relationship.

%Made by Leighton Clarke, June 5 2017

clear all
clc
%Assigned specifications
R1 = 0.45;
R2 = 1.05;
Radj = linspace(0, 10);
X1 = 1.8;
X2 = 2.1;
Xm = 105;
V = 480;
Prated = 15 * 746;
nsync = 1800;
nmlow = 1625;
nmhigh = 1745;

%Plot torque vs resistance at 1625 rpm
s = 1 - nmlow/nsync;
for i = 1:100
    Ztot(i) = R1 + j*X1 + (1/(j*Xm) + 1/((R2+Radj(i))*(1+(1-s)/s) + j*X2))^-1;
    I1(i) = V / Ztot(i);
    I2(i) = I1(i) * j*Xm/((R2+Radj(i))/s + j*(Xm + X2));
    Pdev(i) = 3*I2(i)^2*(R2+Radj(i))*(1-s)/s;
    Tind(i) = Pdev(i)/(nmlow / 60 * 2 * pi);
    P1625(i) = I2(i)^2*Radj(i);
end
plot(Radj,abs(Tind),'k')
hold on

%Plot torque vs resistance at 1745
s = 1 - nmhigh/nsync;
for i = 1:100
    Ztot(i) = R1 + j*X1 + (1/(j*Xm) + 1/((R2+Radj(i))*(1+(1-s)/s) + j*X2))^-1;
    I1(i) = V / Ztot(i);
    I2(i) = I1(i) * j*Xm/((R2+Radj(i))/s + j*(Xm + X2));
    Pdev(i) = 3*I2(i)^2*(R2+Radj(i))*(1-s)/s;
    Tind(i) = Pdev(i)/(nmhigh / 60 * 2 * pi);
    P1745(i) = I2(i)^2*Radj(i);
end
plot(Radj,abs(Tind),'--k')
xlabel('Potentiometer Resistance (ohms)');
ylabel('Induced Torque (N-m)');
title('Torque vs. Potentiometer Resistance');
legend('1625 rpm', '1745 rpm');

%Plot dissipated power vs resistance for both speeds
figure
plot(Radj,abs(P1625),'k')
hold on
plot(Radj,abs(P1745),'--k')
xlabel('Potentiometer Resistance (ohms)');
ylabel('Power Dissipated (W)');
title('Power Dissipated in Potentiometer vs. Resistance Value');
legend('1625 rpm', '1745 rpm')

%Plot four torque speed curves for evenly spaced resistance values
Radj = [0.72,2,3.3,4.6];
nm = linspace(0, 1800);
s = 1 - nm./nsync;
figure
hold on
for k = 1:4
    for i = 1:100
        Ztot(i) = R1 + j*X1 + (1/(j*Xm) + 1/((R2+Radj(k))*(1+(1-s(i))/s(i)) + j*X2))^-1;
        I1(i) = V / Ztot(i);
        I2(i) = I1(i) * j*Xm/((R2+Radj(k))/s(i) + j*(Xm + X2));
        Pdev(i) = 3*I2(i)^2*(R2+Radj(k))*(1-s(i))/s(i);
        Tind(i) = Pdev(i)/(nm(i) / 60 * 2 * pi);
    end
    if k == 1
        plot(nm,abs(Tind),'k')
    elseif k == 2
        plot(nm,abs(Tind),'--k')
    elseif k == 3
        plot(nm,abs(Tind),':k')
    elseif k == 4
        plot(nm,abs(Tind),'-.k')
    end
end

xlabel('Rotor Speed (rpm)');
ylabel('Induced Torque (N-m)');
title('Torque-Speed Curves with Rotor Potentiometer');
legend('0.72 ohms','2.0 ohms','3.3 ohms','4.6 ohms');