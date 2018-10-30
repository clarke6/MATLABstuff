function [Va, Vb, Vc, Ia, Ib, Ic] = RelayFilters(data, CTR, VTR)

%Create arrays for secondary signals
t0 = data(:,1);
va = data(:,2)./VTR;
vb = data(:,3)./VTR;
vc = data(:,4)./VTR;
ia = data(:,5)./CTR;
ib = data(:,6)./CTR;
ic = data(:,7)./CTR;
j = sqrt(-1);
kk=180/pi;

%Calculate RC and preallocate memory for new arrays
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
Va=(1/sqrt(2))*(vax+j*vay).*exp(-j*2*pi*60*t);
Vb=(1/sqrt(2))*(vbx+j*vby).*exp(-j*2*pi*60*t);
Vc=(1/sqrt(2))*(vcx+j*vcy).*exp(-j*2*pi*60*t);
Ia=(1/sqrt(2))*(iax+j*iay).*exp(-j*2*pi*60*t);
Ib=(1/sqrt(2))*(ibx+j*iby).*exp(-j*2*pi*60*t);
Ic=(1/sqrt(2))*(icx+j*icy).*exp(-j*2*pi*60*t);

fprintf('\nDigital Relay Simulation:\n')
fprintf('|Va|: %f V, angle(Va): %f degrees\n', abs(Va(end)), angle(Va(end))*kk)
fprintf('|Vb|: %f V, angle(Vb): %f degrees\n', abs(Vb(end)), angle(Vb(end))*kk)
fprintf('|Vc|: %f V, angle(Vc): %f degrees\n', abs(Vc(end)), angle(Vc(end))*kk)
fprintf('|Ia|: %f A, angle(Ia): %f degrees\n', abs(Ia(end)), angle(Ia(end))*kk)
fprintf('|Ib|: %f A, angle(Ib): %f degrees\n', abs(Ib(end)), angle(Ib(end))*kk)
fprintf('|Ic|: %f A, angle(Ic): %f degrees\n', abs(Ic(end)), angle(Ic(end))*kk)
