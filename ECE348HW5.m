%Leighton Clarke
%ECE 348, HW#5
%May 22, 2017
%This program plots the induced torque, developed power, and line current
%of an induction motor, given a list of specifications. It creates a table
%of motor operating parameters based on the calculated data.

clear all
clc
%Assigned specifications
R1 = 0.45;
R2 = 1.05;
X1 = 1.8;
X2 = 2.1;
Xm = 105;
V = 480;
fe = 60;
Poles = 4;
Prated = 15 * 746;

%Calculate total impedance, currents, developed power, and torque
nsync = 120*fe/Poles;
nm = linspace(0, 1800, 10000);
s = 1 - nm./nsync;
for i = 1:10000
Ztot(i) = R1 + j*X1 + (1/(j*Xm) + 1/(R2*(1+(1-s(i))/s(i)) + j*X2))^-1;
I1(i) = V / Ztot(i);
I2(i) = I1(i) * j*Xm/(R2/s(i) + j*(Xm + X2));
Pdev(i) = 3*I2(i)^2*R2*(1-s(i))/s(i);
Tind(i) = Pdev(i)/(nm(i) / 60 * 2 * pi);
end

%Create plot
plot(nm,abs(Pdev/1000),'k')
hold on
plot(nm,abs(Tind),'--k')
plot(nm,abs(I1),':k')
title('Line current, developed power, and torque vs rotor speed')
xlabel('Rotor speed (rpm)')  
ylabel('Torque/current/power')
legend('Power (kW)','Torque (N-m)','Current (A)','Location','Northwest')

%Parse data for operating parameters
Istart = abs(I1(1));
Inl = abs(I1(10000));
Tpullout = max(abs(Tind));
Pmax = max(abs(Pdev));
Tlockedrotor = abs(Tind(2));

for i = 1:10000
    if (abs(abs(Pdev(i)) - Prated) < 100) && nm(i) > 1000
        nrated = nm(i);
        Irated = abs(I1(i));
        Trated = abs(Tind(i));
        srated = abs(s(i));
    end
end

%Create table
Variable = {'Starting line current';'No-load line current';'Full-load line current';'Locked-rotor torque';'Pullout torque';'Rated torque';'Peak power';'Rated speed';'Slip at rated speed'};
Value = [Istart; Inl; Irated; Tlockedrotor; Tpullout; Trated; Pmax; nrated; srated;];

Unit = categorical({'A';'A';'A';'N-m';'N-m';'N-m';'W';'rpm';'N/A'});
T = table(Value, Unit, 'Rownames', Variable);
disp(T)
