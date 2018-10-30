%{
This program calculates the sending-end voltage of a 300 km, 50 MW,
220 kV transmission line using the long line model.  It does this by
calculating g(amma), Zp(rime), and Yp(rime).
Leighton Clarke
ECE 347, HW 5, text problem 9-16
3-6-2017
%}

clear all
clc

Vrl = 220e3;
Vr = Vrl / sqrt(3);
P = 50e6;
Z = 23 + 75*j;
Y = 500e-6*j;
PF = .88;
d = 300e3;

Ir = (P/PF) / (Vrl * sqrt(3));
IrReal = Ir * cos(-acos(.88));
IrIm = Ir * sin(-acos(.88));
IrCom = IrReal + j*IrIm;
g = sqrt(Z*Y / d^2);
Zp = Z * sinh(g*d) / (g*d);
Yp = Y * tanh(g*d/2) / (g*d/2);

A = Zp*Yp/2 + 1;
B = Zp;

Vs = A*Vr + B*IrCom;
fprintf('Vs = %g + j%g V', real(Vs), imag(Vs))