%MATLAB code for problem 2
%%Part a
clear all
close all
j=sqrt(-1);
kk = 180/pi;
data1 = csvread('fault01_PROB2.csv');
data2 = csvread('fault02_PROB2.csv');
%Values given in problem
CTR = 800/5;
VTR = 500000/sqrt(3)/67;
ZTR = VTR/CTR;
zL0 = 0.1+j*2.6;
zL1 = 0.073+j*0.8;
k0 = (zL0-zL1)/(3*zL1);
ZL1 = zL1*100;
Zr = ZL1*.8/ZTR;
%Use digital relay simulation function to get phasors
Phasors1=RelaySim(data1,CTR,VTR);
%Use Zapp function to get apparent impedances
Zapp1 = Zapp(Phasors1,k0);
%Get Vpol array for cross-polarization
Vpol1 = Vpol(Phasors1);
%Use Zp function to get Zp array
Zp1 = FindZp(Phasors1,Vpol1,k0);
%Get center and radius of MHO characteristic uzing Zr and Zp
Zc1 = zeros(size(Zp1));
r1 = zeros(size(Zp1));
for k = 1:6
    Zc1(:,k) = (Zr+Zp1(:,k))/2;
    r1(:,k) = abs(Zr-Zp1(:,k))/2;
end
%Create arrays for plotting relay characteristics
circle = linspace(0,2*pi);
Xcomp = zeros(length(circle),6);
Ycomp = zeros(length(circle),6);
for k = 1:6
    %Calculate real and imaginary components of circle
    Xcomp(:,k) = r1(end,k)*cos(circle)+real(Zc1(end,k));
    Ycomp(:,k) = r1(end,k)*sin(circle)+imag(Zc1(end,k));
    %Plot circle and apparent impedance for each element
    figure
    plot(Xcomp(:,k),Ycomp(:,k))
    hold on
    grid on
    a=plot(real(Zapp1(:,k)),imag(Zapp1(:,k)),'*');
    plot(real(Zapp1(:,k)),imag(Zapp1(:,k)),'r')
    b=plot([0,real(Zp1(end,k))],[0,imag(Zp1(end,k))],':k');
    c=plot([0,real(Zr)],[0,imag(Zr)],'--k');
    legend([a,c,b],{'Apparent Impedance','Zr','Zp'})
    xlabel('R (secondary)')
    ylabel('X (secondary)')
    switch k
        case 1
            title('Fault 1, A-G')
        case 2
            title('Fault 1, B-G')
        case 3
            title('Fault 1, C-G')
        case 4
            title('Fault 1, A-B')
        case 5
            title('Fault 1, B-C')
        case 6
            title('Fault 1, C-A')
    end
end

%Use MHOout function to calculate MHO relay digital outputs
MHOout1 = MHOout(Phasors1,Vpol1,Zr,k0);

%Plot MHO digital outputs
figure
suptitle('MHO Outputs for Fault 1')
subplot(3,2,1)
plot(Phasors1(:,1),MHOout1(:,1))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOAG')
subplot(3,2,2)
plot(Phasors1(:,1),MHOout1(:,2))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOBG')
subplot(3,2,3)
plot(Phasors1(:,1),MHOout1(:,3))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOCG')
subplot(3,2,4)
plot(Phasors1(:,1),MHOout1(:,4))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOAB')
subplot(3,2,5)
plot(Phasors1(:,1),MHOout1(:,5))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOBC')
xlabel('Time (s)')
subplot(3,2,6)
plot(Phasors1(:,1),MHOout1(:,6))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOCA')
xlabel('Time (s)')

%Repeat all of part a for fault 2
Phasors2=RelaySim(data2,CTR,VTR);
%Use Zapp function to get apparent impedances
Zapp2 = Zapp(Phasors2,k0);
%Get Vpol array for cross-polarization
Vpol2 = Vpol(Phasors2);
%Use Zp function to get Zp array
Zp2 = FindZp(Phasors2,Vpol2,k0);
%Get center and radius of MHO characteristic uzing Zr and Zp
Zc2 = zeros(size(Zp2));
r2 = zeros(size(Zp2));
for k = 1:6
    Zc2(:,k) = (Zr+Zp2(:,k))/2;
    r2(:,k) = abs(Zr-Zp2(:,k))/2;
end
%Create arrays for plotting relay characteristics
circle = linspace(0,2*pi);
Xcomp = zeros(length(circle),6);
Ycomp = zeros(length(circle),6);
for k = 1:6
    %Calculate real and imaginary components of circle
    Xcomp(:,k) = r2(end,k)*cos(circle)+real(Zc2(end,k));
    Ycomp(:,k) = r2(end,k)*sin(circle)+imag(Zc2(end,k));
    %Plot circle and apparent impedance for each element
    figure
    plot(Xcomp(:,k),Ycomp(:,k))
    hold on
    grid on
    a=plot(real(Zapp2(:,k)),imag(Zapp2(:,k)),'*');
    plot(real(Zapp2(:,k)),imag(Zapp2(:,k)),'r');
    b=plot([0,real(Zp2(end,k))],[0,imag(Zp2(end,k))],':k');
    c=plot([0,real(Zr)],[0,imag(Zr)],'--k');
    xlabel('R (secondary)')
    ylabel('X (secondary)')
    legend([a,c,b],{'Apparent Impedance','Zr','Zp'})
    switch k
        case 1
            title('Fault 2, A-G')
        case 2
            title('Fault 2, B-G')
        case 3
            title('Fault 2, C-G')
        case 4
            title('Fault 2, A-B')
        case 5
            title('Fault 2, B-C')
        case 6
            title('Fault 2, C-A')
    end
end

%Use MHOout function to calculate MHO relay digital outputs
MHOout2 = MHOout(Phasors2,Vpol2,Zr,k0);

%Plot MHO digital outputs
figure
suptitle('MHO Outputs for Fault 2')
subplot(3,2,1)
plot(Phasors2(:,1),MHOout2(:,1))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOAG')
subplot(3,2,2)
plot(Phasors2(:,1),MHOout2(:,2))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOBG')
subplot(3,2,3)
plot(Phasors2(:,1),MHOout2(:,3))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOCG')
subplot(3,2,4)
plot(Phasors2(:,1),MHOout2(:,4))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOAB')
subplot(3,2,5)
plot(Phasors2(:,1),MHOout2(:,5))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOBC')
xlabel('Time (s)')
subplot(3,2,6)
plot(Phasors2(:,1),MHOout2(:,6))
ylim([-.25,1.25])
yticks([0,1])
yticklabels({'No op','op'})
ylabel('MHOCA')
xlabel('Time (s)')

%%Part b
%Distance of fault 1
m1 = imag(Zapp1(:,2))*ZTR./imag(ZL1);
figure
plot(Phasors1(:,1),m1)
xlabel('Time (s)')
ylabel('Per-unit distance')
title('m vs. Time for Fault 1')
dist1 = m1(end)*100;
fprintf('\nFault 1 distance: %f miles\n', dist1)

%Calculating distance of fault 2 is redundant since no fault was detected,
%but here is the result of calculating m2 anyways, assuming a C-G fault.
m2 = imag(Zapp2(:,3))*ZTR./imag(ZL1);
figure
plot(Phasors1(:,1),m2)
xlabel('Time (s)')
ylabel('Per-unit distance')
title('m vs. Time for Fault 2')
dist2 = m2(end)*100;
fprintf('Fault 2 distance: %f miles\n', dist2)

%%Part c
%Use provided capacitance and reactance values to calculate traveling wave
%velocity in line
Cmile = 0.013e-6;
Xmile = 0.8;
Lmile = Xmile/(2*pi*60);
Cmeter = Cmile/1609.34;
Lmeter = Lmile/1609.34;
vel = 1/sqrt(Lmeter*Cmeter);
%Calculate time difference between fault 1 wave arrivals at bus 2 and bus 3
dt1 = (dist1*1609.34*2 - 100*1609.34)/vel;
%Calculate theoretical time difference for fault 2 waves at bus 1 and bus
%2, assuming the location is accurate
dt2 = ((50+dist2)*1609.34*2 - 50*1609.34)/vel;
fprintf('\nFault 1 time difference (bus 2 - bus 3): %f seconds', dt1)
fprintf('\nFault 2 time difference (bus 1 - bus 2): %f seconds\n', dt2)

%%
% <include>RelaySim.m</include>
%
%%
% <include>Zapp.m</include>
%
%%
% <include>Vpol.m</include>
%
%%
% <include>FindZp.m</include>
%
%%
% <include>MHOout.m</include>
%