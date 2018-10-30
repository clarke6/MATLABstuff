%Leighton Clarke
%Program to create plot for problem 3

clear all
clc

%Define nominal and inrush current
In = 18.736;
inrush = 18.736 * 12;

%Create arrays for current and time values from provided fuse plots
Im = [1700 510 180 60 41];
Ic = [5000 2000 700 190 70 50];
tm = [0.01 0.1 1 10 100];
tc = [0.02 0.03 0.1 1 10 100];

%Create arrays for currents and times from transformer damage curve
Idpu = [2 3 4 8 40 80 200 400];
Id = Idpu .* In;
td = [950 200 100 20 0.8 0.2 0.03 0.008];

%Plot curves and IEEE point in log-log scale
loglog(Id,td,'k')
hold on
loglog(Im,tm,'--k')
loglog(Ic,tc,':k')
loglog(inrush,0.1,'*')
legend('Transformer Damage Curve', 'Fuse Melting Curve', 'Fuse Total Clearing Curve', 'IEEE Inrush Point')
xlabel('Current (A)')
ylabel('Time (s)')
title('Problem 3: Transformer Damage and Fuse Curves')