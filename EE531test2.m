%This program plots the damage curve of a transformer and the operation
%time curve of a nearby relay. The IEEE inrush point and significant
%currents are calculated and marked.
%Leighton Clarke, Feb. 17 2018

clear all
clc

%Enter defined values and calculate primary currents
Ipui1=35;
TD1=2;
Ipu1=4;
CTR1=400/5;
Ipickup1=Ipu1*CTR1*(34.5/sqrt(3))/230;
Iinst1=Ipui1*CTR1*(34.5/sqrt(3))/230;
Ipui3=50;
TD3=2.7;
Ipu3=4.5;
CTR3=200/5;
tinst=0.025;
In=50E6/(sqrt(3)*230E3);
Iinrush=12*In;
Iinst3=Ipui3*CTR3;
I3pf=5870*(34.5/sqrt(3))/230;
Ipickup3=Ipu3*CTR3;

%Create arrays for currents and times from transformer damage curve
Idpu = [2 3 4 8 40 80 200 400];
Id = Idpu .* In;
td = [950 200 100 20 0.8 0.2 0.03 0.008];

%Create arrays for current axis and operation time values
I = linspace(10,10000,1000);
top1 = [];
top3 = [];

%Use the operation time function for very inverse U3 relay and defined
%instantaneous pickup current to create relay time curve
for k = 1:1000
    if I(k)<=Ipui3*CTR3
        top3(k)=TD3*(0.0963+3.88/((I(k)/(CTR3*Ipu3))^2-1));
    elseif I(k)>Ipui3*CTR3
        top3(k)=tinst;
    end
    
    I1 = I(k) * 230/(34.5/sqrt(3));
    if I1<=Ipui1*CTR1
        top1(k)=TD1*(0.0963+3.88/((I1/(CTR1*Ipu1))^2-1));
    elseif I1>Ipui1*CTR1
        top1(k)=tinst;
    end
end

%Create plot, legend, and labels
loglog(Id,td,'--k')
hold on
loglog(I,top1,'-.k')
loglog(I,top3,'k')
loglog(Iinrush,0.1,'*k')
xticks([10,Ipickup1,100,Ipickup3,Iinst1,I3pf,1000,Iinrush,Iinst3,10000,100000])
xticklabels({10,'Ipu1',1E2,'Ipu3','Ipui1','If',1E3,'Iinrush','Ipui3',1E4,1E5})
legend('Transformer Damage Curve', 'Feeder Relay (1) Time Curve', '230kV Relay (3) Time Curve', 'IEEE Inrush Point')
xtickangle(45)
xlabel('Current at 230kV Relay (A)')
ylabel('Time (s)')
title('Transformer Damage and Relay Coordination')
grid on