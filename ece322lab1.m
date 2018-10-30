x = linspace(-120,0);

    for i = 1:100
        y1(i)=0.0169*x(i) + 1.7846;
        y2(i)=0.0331*x(i) + 3.2468;
        y3(i)=0.0473*x(i) + 4.7074;
        y4(i)=0.0614*x(i) + 6.1674;
        y5(i)=0.0756*x(i) + 7.6061;
        y6(i)=0.0881*x(i) + 8.7397;
       
    end
plot(x,y1)
hold on
plot(x,y2)
hold on
plot(x,y3)
hold on
plot(x,y4)
hold on
plot(x,y5)
hold on
plot(x,y6)
hold on

title('Ic vs. Vce (Extrapolated)')
xlabel('Vce (V)')
ylabel('Ic (mA)')