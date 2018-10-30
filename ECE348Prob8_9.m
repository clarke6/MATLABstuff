%This program plots the torque-speed characteristic for a DC shunt motor
%with and without a 1200 A-turns full-load armature reaction.

%Created by Leighton Clarke, June 5 2017

close all
clear all
clc

%This data comes from the textbook's .m file
If_plot = 0:0.04:2;
Ea_plot = [12.55 27.09 41.36 55.29 68.85 82.00 94.72 106.98 118.76 130.05 140.83 151.10 160.86 170.10 178.84 187.07 194.82 202.08 208.88 215.24 221.16 226.67 231.79 236.54 240.94 245.02 248.80 252.3 255.56 258.58 261.39 264 266.46 268.77 270.95 273.01 274.98 276.86 278.66 280.40 282.08 283.70 285.26 286.77 288.22 289.59 290.88 292.08 293.16 294.1 294.88]';
n0 = 1200;

%Motor parameters
Il = linspace(0,110);
V = 240;
Ra = 0.19;
Rf = 75;
Radj = 175;
Nf = 2700;
f0_Arm = 1200;

%Calculate armature current and induced voltage
Ia = Il * (Radj + Rf) / (Radj + Rf + Ra);
Ea = V - Ia * Ra;

%Calculate armature reaction MMF
Arm = Ia / 110 * f0_Arm;

%Calculate field current with and without arm. reaction
If = V / (Rf + Radj);
If_Arm = If - Arm / Nf;

%Interpolate given data for magnetization curve
Ea0 = interp1(If_plot,Ea_plot,If);
Ea0_Arm = interp1(If_plot,Ea_plot,If_Arm);

%Calculate motor speed with and without arm. reaction
n = (Ea ./ Ea0) * n0;
n_Arm = (Ea ./ Ea0_Arm) * n0;

%Calculate induced with and without arm. reaction
Tind = Ea .* Ia ./ (n * 2 * pi / 60);
Tind_Arm = Ea .* Ia ./ (n_Arm * 2 * pi / 60);

%Plot label curves
plot(Tind, n, 'k');
hold on
plot(Tind_Arm, n_Arm, '--k');
xlabel('Induced Torque (N-m)');
ylabel('Speed (rpm)');
title('Torque-Speed Characteristic of DC Shunt Motor');
legend('No Armature Reaction','With Armature Reaction');