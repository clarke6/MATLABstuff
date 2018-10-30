plot(0,0)
xlabel('x axis')
ylabel('y axis')
hold on
i = 1;

while i <= 5
    x = i;
    y = x^2;
    
   plot(x,y,'s')
   
   i = i + 1
end