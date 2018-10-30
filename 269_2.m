 % Duty Cycle vs. Vg for Cuk Converter

 % Define the plant using an external file
 cuk_parameters;

 % Define the range of input voltage to sweep
 vg = 9:0.01:14;

 d_new = Vo./(vg + Vo);
 ∆_d = d_new−d;

 % Plot change in duty cycle as vg is swept to keep Vo at 24 V
 figure;
 plot(vg, ∆_d);
 xlabel('Input Voltage (V)');
 ylabel('\Delta Duty Cycle');
 axis([9 14 −0.04 0.08]);
 grid on

 % Plot duty cycle as vg is swept to keep Vo at 24 V
 figure;
 plot(vg, d_new);
 xlabel('Input Voltage (V)');
 ylabel('Duty Cycle');
 axis([9 14 0.62 0.74]);
 grid on