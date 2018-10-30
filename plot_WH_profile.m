time = 0:1:23;
volume = [0.5 0 0 0 0 1 3.5 3.5 3.5 3 2.75 2.25 2 1.75 1.5 1.75 2 2.5 3 3 2.75 2.25 1.5 1];
bar(time, volume)
axis([-1 24 0 4])
ax = gca;
ax.XTick = 0:23;
ax.YTick = 0:0.25:4;
grid on
xlabel('Time of day (hours)')
ylabel('Draw Volume (gallons)')
title('Water Heater Daily Usage Profile')