function [ISA, ISB, ISC, IRA, IRB, IRC] = DigiDiffRelay(data, CTR)
    
%Read from data file
t0 = data(:,1);
isa = data(:,2)./CTR;
isb = data(:,3)./CTR;
isc = data(:,4)./CTR;
ira = -data(:,5)./CTR;
irb = -data(:,6)./CTR;
irc = -data(:,7)./CTR;

%Simulate analog RC filter
fc = 600;
RC = 1 / (2*pi*fc);
isaf = zeros(length(t0),1);
isbf = zeros(length(t0),1);
iscf = zeros(length(t0),1);
iraf = zeros(length(t0),1);
irbf = zeros(length(t0),1);
ircf = zeros(length(t0),1);

for k = 2:length(t0)
    isaf(k)=((t0(k)-t0(k-1))*isa(k)+RC*isaf(k-1))/((t0(k)-t0(k-1))+RC);
    isbf(k)=((t0(k)-t0(k-1))*isb(k)+RC*isbf(k-1))/((t0(k)-t0(k-1))+RC);
    iscf(k)=((t0(k)-t0(k-1))*isc(k)+RC*iscf(k-1))/((t0(k)-t0(k-1))+RC);
    iraf(k)=((t0(k)-t0(k-1))*ira(k)+RC*iraf(k-1))/((t0(k)-t0(k-1))+RC);
    irbf(k)=((t0(k)-t0(k-1))*irb(k)+RC*irbf(k-1))/((t0(k)-t0(k-1))+RC);
    ircf(k)=((t0(k)-t0(k-1))*irc(k)+RC*ircf(k-1))/((t0(k)-t0(k-1))+RC);
end

%Interpolate data with lower sample rate
Nsc = 16;
fs = Nsc * 60;
Ts = 1 /fs;
ts = 0:Ts:t0(end);
isafi = interp1(t0, isaf, ts);
isbfi = interp1(t0, isbf, ts);
iscfi = interp1(t0, iscf, ts);
irafi = interp1(t0, iraf, ts);
irbfi = interp1(t0, irbf, ts);
ircfi = interp1(t0, ircf, ts);

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
n = 1:Nsc;
A = (2/Nsc)*cos(2*pi*n/Nsc);
B = (2/Nsc)*sin(2*pi*n/Nsc);
isaRe = filter(A, 1, isafi);
isbRe = filter(A, 1, isbfi);
iscRe = filter(A, 1, iscfi);
iraRe = filter(A, 1, irafi);
irbRe = filter(A, 1, irbfi);
ircRe = filter(A, 1, ircfi);
isaIm = filter(B, 1, isafi);
isbIm = filter(B, 1, isbfi);
iscIm = filter(B, 1, iscfi);
iraIm = filter(B, 1, irafi);
irbIm = filter(B, 1, irbfi);
ircIm = filter(B, 1, ircfi);

%Calculate phasors
ISA = (1/sqrt(2))*(isaRe+1i*isaIm).*exp(-1i*2*pi*60*ts);
ISB = (1/sqrt(2))*(isbRe+1i*isbIm).*exp(-1i*2*pi*60*ts);
ISC = (1/sqrt(2))*(iscRe+1i*iscIm).*exp(-1i*2*pi*60*ts);
IRA = (1/sqrt(2))*(iraRe+1i*iraIm).*exp(-1i*2*pi*60*ts);
IRB = (1/sqrt(2))*(irbRe+1i*irbIm).*exp(-1i*2*pi*60*ts);
IRC = (1/sqrt(2))*(ircRe+1i*ircIm).*exp(-1i*2*pi*60*ts);

end