clc
clear all
close all

%Problem 1
%Create closed-loop transfer function and plot response
%Add one degree to denominator and use "step" to get ramp
n = [-14.11, 345];
d = [1, 8.89, 345, 0];
T1 = tf(n,d);
step(T1)
title('Ramp Response (Problem 1)')

%Problem 2
%Create open loop transfer function
n2 = [80*15, 80*15*2.5];
d2 = [1, 34.28, 292.49, 137.8, 0];
T2 = tf(n2,d2);
%Plot root locus
figure
rlocus(T2)
title('Root Locus (Problem 2)')
%Create closed loop transfer function
figure
%d2cl = [1, 34.28, 292.49, 1668.2, 3826];
%T2cl = tf(n2,d2cl);
T2cl = T2 / (1 + T2);
%Plot step response
step(T2cl)
title('Step Response (Problem 2)')
