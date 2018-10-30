Ra=0.3;
Rf=0.3;
kf=0.025;
Va=100;

w=linspace(0,2000);
T=zeros(1,100);
for i = 1:100
    T(i) = kf*Va^2/(Ra+Rf+kf*w(i)*2*pi/60)^2;
end

plot(w,T,'k')
hold on
ylim([0 800])
plot(w(1),T(1),'ko')
plot(w(90),T(90),'ko')
xlabel('Speed (rpm)')
ylabel('Torque (N-m)')
title('Torque vs. Speed Curve for Series DC Motor')