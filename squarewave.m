clear all
clc
close all
x = linspace(0,20/1);
square = zeros(1,length(x));
for t = 1:100
    for n = 1:2:1001
        square(t)=square(t) + (sin(n*2*pi*60*x(t)))/n;
    end
    square(t) = 4/pi * square(t);
end
plot(x,square)
        