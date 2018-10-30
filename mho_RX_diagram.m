%mho circle characteristics
theta = 0:.01:(2*pi);
% p=0;
% q=-20:0.01:20;
% r=-20:0.01:20;
% s=0;

 
%Zone 1 setting   
 
line_length = 300;
F = 50;
R1 = 0.029792;
L1 = 0.00105678;
zone1_percentage = 80;
zone1_pu = zone1_percentage/100;
angle_Zone1 = atan(2*pi*F*L1/R1);
 
%radius of Zone 1 circle
radius_zone1 = (zone1_pu*line_length*sqrt(R1^2+(2*pi*F*L1)^2))/2;
 
%centre of Zone 1 circle (a,b)
b = (sin(angle_Zone1))*radius_zone1;
a = b/(tan(angle_Zone1));
 
%circle of radius Zone 1 centre at (a,b)
c = radius_zone1*cos(theta)+a;
d = radius_zone1*sin(theta)+b;


%Zone 2 setting 
 
zone2_percentage = 120;
zone2_pu = zone2_percentage/100;
angle_Zone2 = atan(2*pi*F*L1/R1);
 
%radius of Zone 2 circle
radius_zone2 = (zone2_pu*line_length*sqrt(R1^2+(2*pi*F*L1)^2))/2;
 
%centre of Zone 2 circle (e,f)
f = (sin(angle_Zone2))*radius_zone2;
e = f/(tan(angle_Zone2));
 
%circle of radius Zone 2 centre at (e,f)
g = radius_zone2*cos(theta)+e;
h = radius_zone2*sin(theta)+f;

figure(1);
% plot(p,q,r,s,c,d,'b',g,h,'b');
plot(c,d,'k',g,h,'b');
axis([-200 400 -50 400]);
xlabel('Resistance (R) ohm');
ylabel('Reactance (X) ohm');
grid on
hold on
plot(Rsa,Xsa,'r+',R,X,'g+','Markersize',1);
legend('Zone-1','Zone-2','Zc','Zfrc');