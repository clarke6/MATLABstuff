day1 = [0 0 0 0 0 9 25 4 0 0 0 0 0 0 0 0 0 7 7 15 28 0 4 0];
hours = 0:23;
bar(hours, day1)
axis([0 23 0 30])
ax = gca;
ax.XTick = 0:23;
xlabel('Hour of Day')
ylabel('Minutes of WH Activation')
title('WH Usage: Tank 1, 8-1-2017')
