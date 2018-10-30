%Leighton Clarke
%ECE 311, HW #2
%4-15-2017

%Part d
K = (1/1.5)^2;  %Define K for time constant of 1.5
Kv = 2 / sqrt(K);   %Define Kv
NUM = K;    %Numerator of transfer function
DEN = [1, K*Kv, K]; %Denominator of transfer function

T = tf(NUM, DEN);   %Create transfer function
step(T) %Plot step response
title('Figure 1: Step Response')

%Part f
DEN2 = [1, 0, K];   %New denominator for undamped system
T2 = tf(NUM, DEN2); %Create transfer function
figure
time = linspace(0,20);  %Create array for x-axis
step (T2, time) %Plot step response
title('Figure 2: Step Response (Undamped System)')
