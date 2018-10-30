 % SSA Model of Minimum−Phase Cuk Converter

 % Define physical inputs
 clear
 close all

 Vg = 12;
 Vo = 24;
 r = 28;
 c1 = 2*10^−6;
 c2 = 2*10^−5;
 d = Vo/(Vg+Vo);
 l1 = 5*10^−4;
 r1 = 0.01;
 l2 = 7.5*10^−3;
 r2 = 0.01;
 m = −1.5*10^−3;

 % Calculate state variables for each switch position
 s2 = l1*l2−m^2;
 A1 = [−1/(r*c2) 0 1/c2 0;
 0 0 −1/c1 0;
 −l1/s2 l1/s2 −l1*r2/s2 m*r1/s2;
 m/s2 −m/s2 m*r2/s2 −l2*r1/s2];
 B1= [0; 0; −m/s2; l2/s2];
 B2 = B1;
 A2 = [−1/(r*c2) 0 1/c2 0;
 0 0 0 1/c1;
 −l1/s2 m/s2 −l1*r2/s2 m*r1/s2;
 m/s2 −l2/s2 m*r2/s2 −l2*r1/s2];
 C1 = [1 0 0 0];
 C2 = C1;
 D = 0;

 % Calculate combined state variables given duty cycle d
 d1 = 1−d;
 A = d*A1+d1*A2;
 B = d*B1+d1*B2;
 C = d*C1+d1*C2;

 % Calculate B matrix for control input − Bd
 X = −A^−1*B*Vg;
 Bd = (A1−A2)*X;

 % Define the state space model from vg to vo
 SYS = ss(A,B,C,D);

 % Define the state space model from d to vo
 SYSd = ss(A,Bd,C,D);

 % Unit step response from disturbance input vg to vo
 figure;
 [y t] = step(SYS,0.03);
 plot(t, Vo+y);
 xlabel('Time (s)');
 ylabel('Amplitude (V)');
 grid on

 % Pole−zero map from d to vo
 figure;
 pzmap(SYSd);
 h = gcr;
 h.AxesGrid.TitleStyle.FontSize = 10;
 h.AxesGrid.XLabelStyle.FontSize = 10;
 h.AxesGrid.YLabelStyle.FontSize = 10;
 grid on

 % Open−loop Bode Plot from d to vo
 figure;
 margin(SYSd);
 h = gcr;
 h.AxesGrid.Xunits = 'Hz';
 h.AxesGrid.TitleStyle.FontSize = 10;
 h.AxesGrid.XLabelStyle.FontSize = 10;
 h.AxesGrid.YLabelStyle.FontSize = 10;

 clear A1 A2 B1 B2 C1 C2% X s2 r c1 c2 l1 l2 m d1