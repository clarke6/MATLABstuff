clc
clear all
close all

%Problem 3
%Create Bode plot and step response of uncompensated system
nu = [20];
du = [1 7 6];
Tu = tf(nu,du);
margin(Tu)
title('Uncompensated System Bode Diagram (Problem 3)')

T1 = Tu / (1 + Tu);
figure
step(T1)
title('Uncompensated System Step Response (Problem 3)')

%Create Bode plot of compensated system using open-loop transfer function
n = [20  14.2];
d = [1 7 6 0];
T = tf(n,d);
figure
margin(T)
title('Compensated System Bode Diagram (Problem 3)')
%Plot step response using closed-loop transfer function
Tcl = T / (1 + T);
figure 
step(Tcl)
title('Compensated System Step Response (Problem 3)')

%Problem 4
%Create Bode plot using open-loop transfer function
n2 = [140.6 140.6*4.53];
d2 = [1 21.3 106.1 85.8];
T2 = tf(n2,d2);
figure
margin(T2)
title('Compensated System Bode Diagram (Problem 4)')
%Plot step response using closed-loop transfer function
T2cl = T2 / (1 + T2);
figure
step(T2cl)
title('Compensated System Step Response (Problem 4)')