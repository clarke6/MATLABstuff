clear all
close all

%define parameters given in circuit description

Vg = 12;
Vo = 24;
r = 28;
c1 = 2*10^-6;
c2 = 2*10^-5;
d = Vo/(Vg+Vo);
l1 = 5*10^-4;
r1 = 0.01;
l2 = 7.5*10^-3;
r2 = 0.01;
m = -1.5*10^-3;
s2 = l1*l2-m^2;

%Define ss variables for switch on(1) and off(2)

A1=[-1/(r*c2), 0, 1/c2, 0;
    0, 0, -1/c1, 0;
    -l1/s2, l1/s2, -l1*r2/s2, m*r1/s2;
    m/s2, -m/s2, m*r2/s2, -l2*r1/s2];

B1=[0;0;-m/s2;l2/s2];
C1=[1,0,0,0];
D1=0;

A2=[-1/(r*c2), 0, 1/c2, 0;
    0, 0, 0, 1/c1; %This 1/c1 term is negative in book!!!
    -l1/s2, m/s2, -l1*r2/s2, m*r1/s2;
    m/s2, -l2/s2, m*r2/s2, -l2*r1/s2];
B2=B1;
C2=C1;
D2=D1;

%Use duty cycle to create ss variables for averaged model
d1 = 1-d;
A = d*A1+d1*A2;
B = d*B1+d1*B2;
C = d*C1+d1*C2;
D = 0;

%Create Bd variable for control input
X = -A^-1*B*Vg;
Bd = (A1-A2)*X;

%The given converter is a 4th order system, so use ITAE normalized poles for
%4th order

ITAE_4 = [-0.424+1i*1.236, -0.424-1i*1.236, -0.626 + 1i*0.4141, -0.626 - 1i*0.4141];

%Sweep pole scaling multiplier, w0, to find first appropriate value
w0 = [10000:1:10200];
max_error = zeros(size(w0));

for n = 1:length(w0)
    P = w0(n) * ITAE_4;
    %use MATLAB's "place" function to locate poles
    k = place(A,Bd,P);
    %recalculate system with these poles
    Abar = [A-Bd*k];
    sys1 = ss(Abar,B,C,D);
    %determine maximum error for step input
    [y, t] = step(sys1, 0.003);
    max_error(n) = max(abs(y));
end
%{
%plot maximum error as function of frequency multiplier
figure;
h = plot(w0,max_error);
xlabel('Frequency multiplier \omega');
ylabel('Maximum Absolute Error (V)');
grid on
%}
%create array of indexes for all errors <= 0.24 V
z = find(max_error <= 0.24);
%select value of w matching lowest index
w = w0(z(1));

%create system using this scaling factor
P = w * ITAE_4;
k = place(A,Bd,P);
Abar = [A-Bd*k];
sys_fsfb = ss(Abar,B,C,D);
%{
%plot closed-loop system response to step disturbance
[y t x] = step(sys_fsfb, 0.002);
figure()
plot(t, Vo+y);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('System unit step response calculated by MATLAB');
grid on

%Control effort plot shows required change in duty ratio
figure;
plot(t, -k*x'+d);
xlabel('Time (s)');
ylabel('Duty Cycle');
grid on

%Run simulink model
sim('ECE551_proj1_part2');
figure;
plot(simout.time,simout.data);
xlabel('Time (s)');
ylabel('Amplitude (V)');
grid on
title('Output from Simulink simulation');
%}