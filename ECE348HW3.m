fnlA = 61;
fnlB = 61.5;
fnlC = 60.5;
SDA = 0.034;
SDB = 0.03;
SDC = 0.026;

fflA = fnlA/(SDA + 1);
fflB = fnlB/(SDB + 1);
fflC = fnlC/(SDC + 1);

SPA = 3/(fnlA-fflA);
SPB = 3/(fnlB-fflB);
SPC = 3/(fnlC-fflC);

Ptot = 0:0.01:10;
fsys = -(Ptot - (SPA*fnlA+SPB*fnlB+SPC*fnlC)) ./ (SPA+SPB+SPC);

PA = SPA .* (fnlA-fsys);
PB = SPB .* (fnlB-fsys);
PC = SPC .*(fnlC-fsys);

plot(Ptot, PA)
hold on
plot(Ptot, PB)
plot(Ptot, PC)

xlabel('Total Power Supplied to Load (MW)')
ylabel('Individual Generator Contribution (MW)')
title('Distribution of Power Supplied by 3 Generators')
legend('Generator A', 'Generator B', 'Generator C')