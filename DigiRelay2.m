function [ts, Vap, Vbp, Vcp, Iap, Ibp, Icp] = DigiRelay2(data, CTR, VTR)
%This is a function that simulates the analog RC filter, data
%interpolation, and 1-cycle digital filters of a digital relay to produce
%voltage and current phasors from raw data.

%Read from data file
t0 = data(:,1);
Va = data(:,2)./VTR;
Vb = data(:,3)./VTR;
Vc = data(:,4)./VTR;
Ia = data(:,5)./CTR;
Ib = data(:,6)./CTR;
Ic = data(:,7)./CTR;

%Simulate analog RC filter
fc = 600;
RC = 1 / (2*pi*fc);
Vaf = zeros(length(t0),1);
Vbf = zeros(length(t0),1);
Vcf = zeros(length(t0),1);
Iaf = zeros(length(t0),1);
Ibf = zeros(length(t0),1);
Icf = zeros(length(t0),1);

Dt=t0(2)-t0(1);
for k=2:length(t0)
Vaf(k)=(Dt*Va(k)+RC*Vaf(k-1))/(Dt+RC);
Vbf(k)=(Dt*Vb(k)+RC*Vbf(k-1))/(Dt+RC);
Vcf(k)=(Dt*Vc(k)+RC*Vcf(k-1))/(Dt+RC);
Iaf(k)=(Dt*Ia(k)+RC*Iaf(k-1))/(Dt+RC);
Ibf(k)=(Dt*Ib(k)+RC*Ibf(k-1))/(Dt+RC);
Icf(k)=(Dt*Ic(k)+RC*Icf(k-1))/(Dt+RC);
end

%Interpolate data with lower sample rate
Nsc = 16;
fs = Nsc * 60;
Ts = 1 /fs;
ts = 0:Ts:t0(end);
Vafi = interp1(t0, Vaf, ts);
Vbfi = interp1(t0, Vbf, ts);
Vcfi = interp1(t0, Vcf, ts);
Iafi = interp1(t0, Iaf, ts);
Ibfi = interp1(t0, Ibf, ts);
Icfi = interp1(t0, Icf, ts);

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
n = 1:Nsc;
A = (2/Nsc)*cos(2*pi*n/Nsc);
B = (2/Nsc)*sin(2*pi*n/Nsc);
VaRe = filter(A, 1, Vafi);
VbRe = filter(A, 1, Vbfi);
VcRe = filter(A, 1, Vcfi);
IaRe = filter(A, 1, Iafi);
IbRe = filter(A, 1, Ibfi);
IcRe = filter(A, 1, Icfi);
VaIm = filter(B, 1, Vafi);
VbIm = filter(B, 1, Vbfi);
VcIm = filter(B, 1, Vcfi);
IaIm = filter(B, 1, Iafi);
IbIm = filter(B, 1, Ibfi);
IcIm = filter(B, 1, Icfi);

%Calculate phasors
Vap = (1/sqrt(2))*(VaRe+1i*VaIm).*exp(-1i*2*pi*60*ts);
Vbp = (1/sqrt(2))*(VbRe+1i*VbIm).*exp(-1i*2*pi*60*ts);
Vcp = (1/sqrt(2))*(VcRe+1i*VcIm).*exp(-1i*2*pi*60*ts);
Iap = (1/sqrt(2))*(IaRe+1i*IaIm).*exp(-1i*2*pi*60*ts);
Ibp = (1/sqrt(2))*(IbRe+1i*IbIm).*exp(-1i*2*pi*60*ts);
Icp = (1/sqrt(2))*(IcRe+1i*IcIm).*exp(-1i*2*pi*60*ts);

end