x = linspace(0,20/1);
square = zeros(1,length(x))
for t = 1:100
    for n = 1:2:9
        square(t)=sin(n*2*pi*60*x(t));
    end
end
plot(x,square)
        