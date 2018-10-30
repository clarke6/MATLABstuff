clear all
clc
Vrated = 120;   %Primary side voltage
Irated = 0.5;   %Rated current Ip
Rc = 6494;
Req = 25;
Xeq = 26.3;
PF = 0.8;
%Calculate VR for rated load with PF 0.8 lagging
I = 0.4 - 0.3j;
VpFL = 120 + (Req + j*Xeq)*I;
VR = (abs(VpFL)-abs(Vrated)) / abs(VpFL) * 100;
%Plot efficiency for power sweep from 0-120%
x = linspace(0,72*PF);
eff = zeros(100);
for t = 1:100
    Pcu = abs((pol2cart(-acos(PF),x(t)/PF)/Vrated))^2 * Req;
    VpFL = 120 + (Req + j*Xeq)*(pol2cart(-acos(PF),x(t)/PF)/Vrated);
    Pcore = abs(VpFL)^2 / Rc;
    Ploss = Pcu + Pcore;
    eff(t) = x(t) / (x(t) + Ploss) * 100;
end

plot(x,eff);
title('8341 Transformer Efficiency vs Load');
xlabel('Load Power (W)');
ylabel('Efficiency (%)');
fprintf('With full load at 0.8 lagging PF, VR = %.1f%%\n', VR);