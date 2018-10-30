function M = RelaySim(data, CTR, VTR)

j = sqrt(-1);
%Create array for outputs
RCfil=zeros(size(data));
data(:,2:4) = data(:,2:4)/VTR;
data(:,5:7) = data(:,5:7)/CTR;

%Calculate RC
fc = 600;
RC = 1 / (2*pi*fc);

%Use recursive equation to simulate analog low-pass filter
for n = 2:length(RCfil(1,:))
    for k = 2:length(RCfil(:,1))
        RCfil(k,n) = ((data(k,1)-data(k-1,1))*data(k,n)+RC*RCfil(k-1,n))/((data(k,1)-data(k-1,1))+RC);
    end
end

%Calculate sampling period and interpolate data
Nsc=16;
fs=Nsc*60;
Dt=1/fs;
t=0:Dt:data(end,1);
Sam = zeros(length(t),7);
Sam(:,1)=t;
for k = 2:7
    Sam(:,k) = interp1(data(:,1), RCfil(:,k),t);
end

%Use sine and cosine digital filters to separate real and imaginary parts
%of fundamental frequency signal
k=1:Nsc;
bc=(2/Nsc)*cos(2*pi*k/Nsc);
bs=(2/Nsc)*sin(2*pi*k/Nsc);
digfilx = zeros(length(t),7);
digfily = zeros(length(t),7);
for k = 2:7
    digfilx(:,k) = filter(bc,1,Sam(:,k));
    digfily(:,k) = filter(bs,1,Sam(:,k));
end

%Create phasor array for each signal
M = zeros(length(t),length(data(1,:)));
M(:,1)=t;
for k = 2:length(M(1,:))
    M(:,k)=(1/sqrt(2))*(digfilx(:,k)+j*digfily(:,k)).*exp(-j*2*pi*60*M(:,1));
end
