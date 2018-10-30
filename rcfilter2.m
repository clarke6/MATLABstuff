function xf=rcfilter2(fc,x,t)
% Home-made RC filter function
% x = input vector
% t = time vector
% fc =cut-off frequence in Hz

RC=1/(2*pi*fc);
Dt=t(2)-t(1);
xf(1:length(t),1)=0;
for k=2:length(t)
xf(k)=(Dt*x(k)+RC*xf(k-1))/(Dt+RC);
end

end