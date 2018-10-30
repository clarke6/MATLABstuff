ECE551_Project_part1;
sim('ECE551_proj1')
figure;
plot(simout.time, simout.data);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Step response of open-loop model in Simulink');
grid on;