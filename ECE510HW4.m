close all
x = [0 0 45 45 90 90 180 180 225 225 270 270 360 360 405];
Nr = [0 1 1 2 2 3 3 2 2 1 1 0 0 1 1];
plot(x,Nr)
xticks([0 45 90 135 180 225 270 315 360]);
yticks([1 2 3]);
yticklabels({'Nr/3','2Nr/3','Nr'});
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ylim([-.5,4]);
xlim([-40,380]);
title('Rotor Turns Function');
ylabel('n_r(\theta)')
xlabel('\theta')

figure(2)
x = [0 45 45 225 225 360];
Ns = [0 0 1 1 0 0];
plot(x,Ns);
xticks([0 45 90 135 180 225 270 315 360]);
yticks([1]);
yticklabels({'Ns'});
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ylim([-1/8, 1.25]);
xlim([-40,380]);
title('Stator Turns Function');
ylabel('n_s(\theta)')
xlabel('\theta')

figure(3)
x = [0 0 45 45 90 90 180 180 225 225 270 270 360 360 405];
Nr = Nr - 1.5;
plot(x,Nr)
xticks([0 45 90 135 180 225 270 315 360]);
yticks([-1.5,-1/2,1/2,1.5]);
yticklabels({'-Nr/2','-Nr/6','Nr/6','Nr/2'});
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ylim([-2,2]);
xlim([-40,370]);
title('Rotor Windings Function');
ylabel('N_r(\theta)')
xlabel('\theta')

figure(4)
x = [0 45 45 225 225 405];
Ns = Ns - 1/2;
plot(x,Ns);
xticks([0 45 90 135 180 225 270 315 360]);
yticks([-1/2,1/2]);
yticklabels({'-Ns/2','Ns/2'});
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ylim([-1, 1]);
xlim([-40,370]);
title('Stator Windings Function');
ylabel('N_s(\theta)')
xlabel('\theta')

figure(5)
subplot(3,1,1);
nA = [0 6 6 0 0 6 6 0 0 6 6 0 0];
x = [0 0 30 30 60 60 90 90 120 120 150 150 180];
plot(x,nA)
title('A-Phase Turns Function');
ylabel('n_A(\theta)')
xlabel('\theta')
ylim([-2,8]);
xticks([0 30 60 90 120 150 180]);
yticks([0 6])

subplot(3,1,2);
nB = [0 0 -6 -6 0 0 -6 -6 0 0 -6 -6 0 0];
x = [0 10 10 40 40 70 70 100 100 130 130 160 160 190];
plot(x,nB)
title('B-Phase Turns Function');
ylabel('n_B(\theta)')
xlabel('\theta')
ylim([-8,2]);
xticks([0 10 40 70 100 130 160 180]);
xlim([0,180]);
yticks([-6,0])

subplot(3,1,3);
nC = [0 0 6 6 0 0 6 6 0 0 6 6 0 0];
x = [0 20 20 50 50 80 80 110 110 140 140 170 170 200];
plot(x,nC)
title('C-Phase Turns Function');
ylabel('n_C(\theta)')
xlabel('\theta')
ylim([-2,8]);
xticks([0 20 50 80 110 140 170 180]);
xlim([0,180]);
yticks([0,6])

figure(6)
subplot(3,1,1);
nA = [0 6 6 0 0 6 6 0 0 6 6 0 0]-3;
x = [0 0 30 30 60 60 90 90 120 120 150 150 180];
plot(x,nA)
title('A-Phase Winding Function');
ylabel('N_A(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 30 60 90 120 150 180]);
yticks([-3,3]);

subplot(3,1,2);
nB = [0 0 -6 -6 0 0 -6 -6 0 0 -6 -6 0 0]+3;
x = [0 10 10 40 40 70 70 100 100 130 130 160 160 190];
plot(x,nB)
title('B-Phase Winding Function');
ylabel('N_B(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 10 40 70 100 130 160 180]);
xlim([0,180]);
yticks([-3,3]);

subplot(3,1,3)
nC = [0 0 6 6 0 0 6 6 0 0 6 6 0 0]-3;
x = [0 20 20 50 50 80 80 110 110 140 140 170 170 200];
plot(x,nC)
title('C-Phase Winding Function');
ylabel('N_C(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 20 50 80 110 140 170 180]);
xlim([0,180]);
yticks([-3,3]);

figure(7)
subplot(3,1,1);
nA = [0 6 6 0 0 6 6 0 0 6 6 0 0]-3;
x = [0 0 30 30 60 60 90 90 120 120 150 150 180];
plot(x,nA,'--')
hold on
deg = linspace(0,180);
plot(deg,4/pi*3*sin(6*deg*pi/180))
title('A-Phase Winding Function vs fundamental component');
ylabel('N_A(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 30 60 90 120 150 180]);
yticks([-3,3]);

subplot(3,1,2);
nB = [0 0 -6 -6 0 0 -6 -6 0 0 -6 -6 0 0]+3;
x = [0 10 10 40 40 70 70 100 100 130 130 160 160 190];
plot(x,nB,'--')
hold on
plot(deg,4/pi*3*sin(6*((deg+20)*pi/180)))
title('B-Phase Winding Function vs fundamental component');
ylabel('N_B(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 10 40 70 100 130 160 180]);
xlim([0,180]);
yticks([-3,3]);

subplot(3,1,3)
nC = [0 0 6 6 0 0 6 6 0 0 6 6 0 0]-3;
x = [0 20 20 50 50 80 80 110 110 140 140 170 170 200];
plot(x,nC,'--')
hold on
plot(deg,4/pi*3*sin(6*((deg-20)*pi/180)))
title('C-Phase Winding Function vs fundamental component');
ylabel('N_C(\theta)')
xlabel('\theta')
ylim([-4,4]);
xticks([0 20 50 80 110 140 170 180]);
xlim([0,180]);
yticks([-3,3]);