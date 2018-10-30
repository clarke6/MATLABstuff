%ECE331, Homework 1
%Leighton Clarke, 4-10-17
%This program finds the closed-loop transfer function of a
%feedback system, plots the system's step response, and
%shows the final value of the step response.

clear all
clc
fprintf('Part a: closed-loop transfer function')
N1 = [1];
N2 = [1 2];
D1 = [1 1];
D2 = [1 3];
M1 = tf(N1,D1); %Transfer function for first element
M2 = tf(N2,D2); %Transfer function for second element
M = series(M1, M2); %Series combination of elements M1 and M2
Transfer = feedback(M, 1)  %Closed-loop transfer function of M1 and M2

step(Transfer)  %Plot transfer function
t = linspace(0,10); %Create time array
step = step(Transfer, t);   %Calculate step response for 0-10
fprintf('Final Value: %g\n', step(100)) %Show value of step response at 10