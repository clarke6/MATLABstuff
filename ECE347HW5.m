%{
This program calculates the per unit distance and net impedances
for a 150 MVA, 500 kV, 95 mile transmission line using 3 suitable
ACSR conductors.  It also calculates the ABCD matrix coefficients.
Leighton Clarke
ECE 347, HW 5, Problem 1
March 1, 2017
%}

clear all
clc

%Provided data
LineDist = 95 * 5.28;   %kilofeet
GMR = [0.0101, 0.0113, 0.0127] ; %feet
GMD = 8;  %feet
r = [0.259, .206, .163];  %ohms/kilofoot
xl = [0.1100, 0.1068, 0.1040]; %ohms/kft
xc = [0.6785, 0.66, 0.6421]; %Megaohms/kft

for i = 1:3
    % Calculate impedances
    Xc = xc(i) / LineDist * log(GMD/GMR(i)); %Megaohms
    yc = j / xc(i); %uS/kft
    Yc = j / Xc;    %uS
    R = r(i) * LineDist;    %Ohms
    Xl = xl(i) * LineDist * log(GMD/GMR(i));  %Ohms
    Z = R + Xl*j;

    % Calculate ABCD matrix coefficients
    A = 1 + Z*Yc/2;
    B = Z;
    C = Yc*(1 + Z*Yc/4);
    D = A;
    
    if i == 1;
        fprintf('Sparrow:\n');
    elseif i == 2;
        fprintf('Robin:\n');
    elseif i == 3;
        fprintf('Raven:\n');
    end
    fprintf('Impedances\n')
    fprintf('r = %.3g ohm/kft\nR = %.3g ohm\nxl = %.3g ohm/kft\n', r(i), R, xl(i))
    fprintf('Xl = %.3g ohm\ny = j%.3g uS/kft\nY = j%.3g uS\n', Xl, imag(yc), imag(Yc))
    fprintf('\nMatrix Coefficients\n')
    fprintf('A = %.3g + j%.3g\nB = %.3g + j%.3g\n', real(A), imag(A), R, Xl)
    fprintf('C = %.3g + j%.3g\nD = %.3g + j%.3g\n\n', real(C), imag(C), real(D), imag(D))
end