figure
x = [75, 125, 200];
y = [7.06, 10.68, 14.26];
plot(x,y,'-*')
xlabel('BIL (kV)')
ylabel('Bushing Height (inches)')
title('Bushing Height vs. BIL for EX-7L Capacitors')
txt1 = '(75, 7.06)';
text(x(1),y(1),txt1)
txt2 = '(125, 10.68)';
txt3 = '(200, 14.26)';
text(x(2),y(2),txt2)
text(x(3),y(3),txt3)

