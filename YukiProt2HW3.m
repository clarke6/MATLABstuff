%Yuki Kato
%Protections 2, homework 3


clear all
close all


%Read data and create arrays for secondary signals
data = csvread('prot2_hw2_2108.csv');
PTR = 525000/115;
CTR = 800/5;
kk=180/pi;
t0 = data(:,1);
va = data(:,2)./PTR;
vb = data(:,3)./PTR;
vc = data(:,4)./PTR;
ia = data(:,5)./CTR;
ib = data(:,6)./CTR;
ic = data(:,7)./CTR;
j = sqrt(-1);
theta = 0:0.01:(2*pi);
VTR = 525000/115;
CTR = 800/5;
ZTR = VTR/CTR;

%Calculate RC and preallocate memory for new arrays
fc=600;
fc = 600;
RC = 1 / (2*pi*fc);
vaf = zeros(length(t0),1);
vbf = zeros(length(t0),1);
vcf = zeros(length(t0),1);
iaf = zeros(length(t0),1);
ibf = zeros(length(t0),1);
icf = zeros(length(t0),1);
%Use recursive equation to simulate analog low-pass filter
for k = 2:length(t0)
    vaf(k)=((t0(k)-t0(k-1))*va(k)+RC*vaf(k-1))/((t0(k)-t0(k-1))+RC);
    vbf(k)=((t0(k)-t0(k-1))*vb(k)+RC*vbf(k-1))/((t0(k)-t0(k-1))+RC);
    vcf(k)=((t0(k)-t0(k-1))*vc(k)+RC*vcf(k-1))/((t0(k)-t0(k-1))+RC);
    iaf(k)=((t0(k)-t0(k-1))*ia(k)+RC*iaf(k-1))/((t0(k)-t0(k-1))+RC);
    ibf(k)=((t0(k)-t0(k-1))*ib(k)+RC*ibf(k-1))/((t0(k)-t0(k-1))+RC);
    icf(k)=((t0(k)-t0(k-1))*ic(k)+RC*icf(k-1))/((t0(k)-t0(k-1))+RC);
end

%Calculate sampling period and interpolate data
Nsc=16;
fs=Nsc*60;
Dt=1/fs;
t=[0:Dt:t0(length(t0))];
vafr=interp1(t0,vaf,t);
vbfr=interp1(t0,vbf,t);
vcfr=interp1(t0,vcf,t);
iafr=interp1(t0,iaf,t);
ibfr=interp1(t0,ibf,t);
icfr=interp1(t0,icf,t);

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
k=[1:Nsc];
bc=(2/Nsc)*cos(2*pi*k/Nsc);
bs=(2/Nsc)*sin(2*pi*k/Nsc);
vax=filter(bc,1,vafr);
vbx=filter(bc,1,vbfr);
vcx=filter(bc,1,vcfr);
vay=filter(bs,1,vafr);
vby=filter(bs,1,vbfr);
vcy=filter(bs,1,vcfr);
iax=filter(bc,1,iafr);
ibx=filter(bc,1,ibfr);
icx=filter(bc,1,icfr);
iay=filter(bs,1,iafr);
iby=filter(bs,1,ibfr);
icy=filter(bs,1,icfr);

%Create phasor array for each signal
j=sqrt(-1);
Va=(1/sqrt(2))*(vax+j*vay).*exp(-j*2*pi*60*t);
Vb=(1/sqrt(2))*(vbx+j*vby).*exp(-j*2*pi*60*t);
Vc=(1/sqrt(2))*(vcx+j*vcy).*exp(-j*2*pi*60*t);
Ia=(1/sqrt(2))*(iax+j*iay).*exp(-j*2*pi*60*t);
Ib=(1/sqrt(2))*(ibx+j*iby).*exp(-j*2*pi*60*t);
Ic=(1/sqrt(2))*(icx+j*icy).*exp(-j*2*pi*60*t);

%Calculate phasors using traditional method
E=525000/sqrt(3);
Zs1=2+30*j;
Zs0=Zs1;
zL1=0.073+0.8*j;
zL0=0.1+2.3*j;
I1=E/(Zs1+zL1*20+Zs1+zL1*20+Zs0+zL0*20);
I2=I1;I0=I1;
V1=E-Zs1*I1;
V2=-Zs1*I2;
V0=-Zs0*I0;
a=1*exp(j*120/kk);
Iar=(I0+I1+I2)/CTR;
Ibr=(I0+(a^2)*I1+a*I2)/CTR;
Icr=(I0+a*I1+(a^2)*I2)/CTR;
Var=(V0+V1+V2)/PTR;
Vbr=(V0+(a^2)*V1+a*V2)/PTR;
Vcr=(V0+a*V1+(a^2)*V2)/PTR;
Ires=Iar+Ibr+Icr;

%MHO relay characteristic
ZL1 = (0.073+j*0.8)*100;
Zr = ZL1 * 0.8;
ZL1s = ZL1/ZTR;
Zrs = Zr/ZTR;
r = abs(Zrs)/2;
rR = r*cos(theta) + real(Zrs)/2;
rX = r*sin(theta) + imag(Zrs)/2;

%Calculate impedances
Zab = (Var-Vbr)/(Iar-Ibr);
Zbc = (Vbr-Vcr)/(Ibr-Icr);
Zca = (Vcr-Var)/(Icr-Iar);
k0 = (zL0-zL1)/(3*zL1);
Zag = Var/(Iar+k0*Ires);
Zbg = Vbr/(Ibr+k0*Ires);
Zcg = Vcr/(Icr+k0*Ires);

%Print results
fprintf('\nDeliverables #2: traditionally calculated values:\n')
fprintf('Va: %f V, angle %f deg\n', abs(Var), angle(Var)*kk)
fprintf('Vb: %f V, angle %f deg\n', abs(Vbr), angle(Vbr)*kk)
fprintf('Vc: %f V, angle %f deg\n', abs(Vcr), angle(Vcr)*kk)
fprintf('Ia: %f A, angle %f deg\n', abs(Iar), angle(Iar)*kk)
fprintf('Ib: %f A, angle %f deg\n', abs(Ibr), angle(Ibr)*kk)
fprintf('Ic: %f A, angle %f deg\n', abs(Icr), angle(Icr)*kk)
fprintf('Zab: %f + j%f secondary ohms\n', real(Zab), imag(Zab))
fprintf('Zbc is undefined\n')
fprintf('Zca: %f + j%f secondary ohms\n', real(Zca), imag(Zca))
fprintf('Zag: %f + j%f secondary ohms\n', real(Zag), imag(Zag))
fprintf('Zbg: %f + j%f secondary ohms\n', real(Zbg), imag(Zbg))
fprintf('Zcg: %f + j%f secondary ohms\n', real(Zcg), imag(Zcg))

%Plot faults and protective zone
figure(1)
plot(rR,rX)
grid on
hold on
plot(real(Zag),imag(Zag),'xk')
plot([0,real(ZL1s)],[0,imag(ZL1s)],'--')
legend('Relay Characteristic','Zag','ZL1')
legend('Location','northwest')
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Deliverables #3: AG Relay (self-polarized, calculated)')

%Create arrays of apparent impedances as seen by relay
Zabr = (Va-Vb)./(Ia-Ib);
Zbcr = (Vb-Vc)./(Ib-Ic);
Zcar = (Vc-Va)./(Ic-Ia);
Iresr = Ia+Ib+Ic;
Zagr = Va./(Ia+k0*Iresr);
Zbgr = Vb./(Ib+k0*Iresr);
Zcgr = Vc./(Ic+k0*Iresr);

%Print results
fprintf('\nDeliverables #4: values measured by digital relay:\n')
fprintf('Va: %f V, angle %f deg\n', abs(Va(end)), angle(Va(end))*kk)
fprintf('Vb: %f V, angle %f deg\n', abs(Vb(end)), angle(Vb(end))*kk)
fprintf('Vc: %f V, angle %f deg\n', abs(Vc(end)), angle(Vc(end))*kk)
fprintf('Ia: %f A, angle %f deg\n', abs(Ia(end)), angle(Ia(end))*kk)
fprintf('Ib: %f A, angle %f deg\n', abs(Ib(end)), angle(Ib(end))*kk)
fprintf('Ic: %f A, angle %f deg\n', abs(Ic(end)), angle(Ic(end))*kk)
fprintf('Zab: %f + j%f secondary ohms\n', real(Zabr(end)), imag(Zabr(end)))
fprintf('Zbc: %f + j%f secondary ohms\n', real(Zbcr(end)), imag(Zbcr(end)))
fprintf('Zca: %f + j%f secondary ohms\n', real(Zcar(end)), imag(Zcar(end)))
fprintf('Zag: %f + j%f secondary ohms\n', real(Zagr(end)), imag(Zagr(end)))
fprintf('Zbg: %f + j%f secondary ohms\n', real(Zbgr(end)), imag(Zbgr(end)))
fprintf('Zcg: %f + j%f secondary ohms\n', real(Zcgr(end)), imag(Zcgr(end)))

%Plot evolution of apparent impedance Zag in complex plane
figure(2)
plot(real(Zagr),imag(Zagr))
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Deliverables #5: Evolution of Zag')

%Plot steady-state of Zag in MHO characteristic
figure(3)
plot(rR,rX)
hold on
plot(real(Zagr(end)),imag(Zagr(end)),'xk')
plot([0,real(ZL1s)],[0,imag(ZL1s)],'--')
grid on
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Deliverables #6: AG Relay (self-polarized, measured)')
legend('Relay Characteristic','Zag','ZL1')
legend('Location','northwest')

%Create digital MHO variables for each datapoint
Sab = (Ia-Ib)*Zrs - (Va-Vb);
Sbc = (Ib-Ic)*Zrs - (Vb-Vc);
Sca = (Ic-Ia)*Zrs - (Vc-Va);
Sag = (Ia+k0*Iresr)*Zrs - Va;
Sbg = (Ib+k0*Iresr)*Zrs - Vb;
Scg = (Ic+k0*Iresr)*Zrs - Vc;

MHOAB = zeros(1,length(Sab));
MHOBC = zeros(1,length(Sab));
MHOCA = zeros(1,length(Sab));
MHOAG = zeros(1,length(Sab));
MHOBG = zeros(1,length(Sab));
MHOCG = zeros(1,length(Sab));

for m = 1:length(MHOAB)
    if real(Sab(m)*conj(Va(m)-Vb(m))) > 0
        MHOAB(m) = 1;
    end
    if real(Sbc(m)*conj(Vb(m)-Vc(m))) > 0
        MHOBC(m) = 1;
    end
    if real(Sca(m)*conj(Vc(m)-Va(m))) > 0
        MHOCA(m) = 1;
    end
    if real(Sag(m)*conj(Va(m))) > 0
        MHOAG(m) = 1;
    end
    if real(Sbg(m)*conj(Vb(m))) > 0
        MHOBG(m) = 1;
    end
    if real(Scg(m)*conj(Vc(m))) > 0
        MHOCG(m) = 1;
    end
end

%Plot digital operation flags
figure('position',[200 200 700 700])
subplot(5,2,[1,2])
plot(t0,va)
title('Deliverables #7: digital MHO operation variables (self-polarized)')
hold on
plot(t0,vb)
plot(t0,vc)
ylabel('V')
subplot(5,2,[3,4])
plot(t0,ia)
hold on
plot(t0,ib)
plot(t0,ic)
ylabel('A')
subplot(5,2,5)
plot(t,MHOAB)
ylabel('MHOAB')
subplot(5,2,6)
plot(t,MHOBC)
ylabel('MHOBC')
subplot(5,2,7)
plot(t,MHOCA)
ylabel('MHOCA')
subplot(5,2,8)
plot(t,MHOAG)
ylabel('MHOAG')
subplot(5,2,9)
plot(t,MHOBG)
ylabel('MHOBG')
subplot(5,2,10)
plot(t,MHOCG)
ylabel('MHOCG')

%Calculate cross-polarized MHO characteristic for traditionally calculated
%values (AG relay only)
Vpolr = (Vbr-Vcr)*exp(j*90/kk);
Vpr = Var - Vpolr;
Zpr = Vpr/(Iar+k0*Ires);
C = (Zrs + Zpr)/2;
rcross = abs(Zrs-Zpr)/2;
rRcross = rcross*cos(theta)+real(C);
rXcross = rcross*sin(theta)+imag(C);

%Plot both characteristics and apparent impedance Zag
figure(5)
plot(rRcross,rXcross)
hold on
plot(rR,rX,'--')
plot(real(Zag),imag(Zag),'xk')
plot([real(Zpr),real(Zrs)],[imag(Zpr),imag(Zrs)])
legend('Cross-polarized characteristic','Self-polarized characteristic','Zag','Zr-Zp')
xlim([-3,3])
ylim([-2,4])
grid on
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Deliverables #8: AG Relay (self and cross polarized, calculated values)')

%Calculate MHO operation variables using cross-polarizing scheme
MHOABcross = zeros(1,length(Sab));
MHOBCcross = zeros(1,length(Sab));
MHOCAcross = zeros(1,length(Sab));
MHOAGcross = zeros(1,length(Sab));
MHOBGcross = zeros(1,length(Sab));
MHOCGcross = zeros(1,length(Sab));

for m = 1:length(MHOABcross)
    if real(Sab(m)*conj(Vc(m)*exp(-j*90/kk))) > 0
        MHOABcross(m) = 1;
    end
    if real(Sbc(m)*conj(Va*exp(-j*90/kk))) > 0
        MHOBCcross(m) = 1;
    end
    if real(Sca(m)*conj(Vb(m)*exp(-j*90/kk))) > 0
        MHOCAcross(m) = 1;
    end
    if real(Sag(m)*conj((Vb(m)-Vc(m))*exp(j*90/kk))) > 0
        MHOAGcross(m) = 1;
    end
    if real(Sbg(m)*conj((Vc(m)-Va(m))*exp(j*90/kk))) > 0
        MHOBGcross(m) = 1;
    end
    if real(Scg(m)*conj((Va(m)-Vb(m))*exp(j*90/kk))) > 0
        MHOCGcross(m) = 1;
    end
end

%Plot MHO variables for cross-polarizing scheme
figure('position',[200 200 700 700])
subplot(5,2,[1,2])
plot(t0,va)
title('Deliverables #9: digital MHO operation variables (cross-polarized)')
hold on
plot(t0,vb)
plot(t0,vc)
ylabel('V')
subplot(5,2,[3,4])
hold on
plot(t0,ia)
plot(t0,ib)
plot(t0,ic)
ylabel('A')
subplot(5,2,5)
plot(t,MHOABcross)
ylabel('MHOAB')
subplot(5,2,6)
plot(t,MHOBCcross)
ylabel('MHOBC')
subplot(5,2,7)
plot(t,MHOCAcross)
ylabel('MHOCA')
subplot(5,2,8)
plot(t,MHOAGcross)
ylabel('MHOAG')
subplot(5,2,9)
plot(t,MHOBGcross)
ylabel('MHOBG')
subplot(5,2,10)
plot(t,MHOCGcross)
ylabel('MHOCG')

%Plot AG element MHO characteristics and final value of Zag at relay
Vpolp = (Vb-Vc).*exp(j*90/kk);
Vpp = Va - Vpolp;
Zpp = Vpp ./ (Ia+k0*Iresr);
Cp = (Zrs + Zpp)/2;
rcrossp = abs(Zrs-Zpp)/2;
rRcrossp = rcrossp(end)*cos(theta)+real(Cp(end));
rXcrossp = rcrossp(end)*sin(theta)+imag(Cp(end));
figure(7)
plot(rRcrossp,rXcrossp)
hold on
plot(rR,rX,'--')
plot(real(Zagr(end)),imag(Zagr(end)),'kx')
plot([real(Zpp(end)),real(Zrs)],[imag(Zpp(end)),imag(Zrs)])
grid on
xlim([-3,3])
ylim([-2,4])
xlabel('R (secondary ohms)')
ylabel('X (secondary ohms)')
title('Deliverables #10: AG relay (self and cross polarized, measured values)')
legend('Cross-polarized characteristic','Self-polarized characteristic','Zag','Zr-Zp')
