%% 
% Advanced Protections, Homework 1
% 
% Leighton Clarke
% 
% April 5, 2018
% 
% 
% 
% Part b (part a not included)
% 
% Read the data file (must be headless and named EE537HW1.csv):

data = csvread('EE537HW1.csv');
%% 
% Extract each column to define the time, voltage, and current arrays:

t = data(:,1);
Va = data(:,2);
Vb = data(:,3);
Vc = data(:,4);
Ia = data(:,5);
Ib = data(:,6);
Ic = data(:,7);
%% 
% Plot voltages vs. time:

figure
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Voltages')
hold on
plot(t,Va)
plot(t,Vb)
plot(t,Vc)
legend('Va', 'Vb', 'Vc')
%% 
% Plot currents vs. time:

figure
xlabel('Time (s)')
ylabel('Current (A)')
title('Currents')
hold on
plot(t,Ia)
plot(t,Ib)
plot(t,Ic)
legend('Ia', 'Ib', 'Ic')
%% 
% Part c
% 
% Calculate the RC value for a simple analog RC filter with a cutoff frequency 
% of 500Hz:
% 
% $$\mathrm{RC}=\backslash \mathrm{frac}\left\lbrace 1\right\rbrace \left\lbrace 
% 2undefined\backslash \mathrm{pi}undefined\mathrm{fc}\right\rbrace$$