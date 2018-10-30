function bcfault(Es,Zs1,Zs0,Er,Zr1,Zr0,ZL1,ZL0,m,Rf,CTR,VTR)
kk=180/pi;
a=1*exp(1i*120/kk);
Zx1=m*ZL1+Zs1;
Zy1=(1-m)*ZL1+Zr1;
Zx0=m*ZL0+Zs0;
Zy0=(1-m)*ZL0+Zr0;
Z1 = Zx1*Zy1/(Zx1+Zy1);
Z0 = Zx0*Zy0/(Zx0+Zy0);
C1 = Zy1/(Zx1+Zy1);
C0 = Zy0/(Zx0+Zy0);
Eth = C1*Es+(1-C1)*Er;
Ipf = (Es-Er)/(Zs1+ZL1+Zr1);
If1 = Eth/(2*Z1+Rf);
alpha = Ipf/If1;
I1 = If1*(C1+alpha);
I2 = -C1*If1;
Ia = (I1+I2)/CTR;
Ib = (a^2*I1 + a*I2)/CTR;
Ic = (a*I1 + a^2*I2)/CTR;
V1 = If1*((C1+alpha)*m*ZL1+Z1+Rf);
V2 = If1*(-C1*m*ZL1+Z1);
Va = (V1+V2)/VTR;
Vb = (a^2*V1 + a*V2)/VTR;
Vc = (a*V1 + a^2*V2)/VTR;

fprintf('\nTraditional Calculations:\n')
fprintf('|Va|: %f V, angle(Va): %f degrees\n', abs(Va(end)), angle(Va(end))*kk)
fprintf('|Vb|: %f V, angle(Vb): %f degrees\n', abs(Vb(end)), angle(Vb(end))*kk)
fprintf('|Vc|: %f V, angle(Vc): %f degrees\n', abs(Vc(end)), angle(Vc(end))*kk)
fprintf('|Ia|: %f A, angle(Ia): %f degrees\n', abs(Ia(end)), angle(Ia(end))*kk)
fprintf('|Ib|: %f A, angle(Ib): %f degrees\n', abs(Ib(end)), angle(Ib(end))*kk)
fprintf('|Ic|: %f A, angle(Ic): %f degrees\n', abs(Ic(end)), angle(Ic(end))*kk)