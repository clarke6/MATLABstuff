ECE551_Project_part2;
sim('ECE551_proj1_part2');
figure;
plot(simout.time, simout.data);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Step response of FSFB model in Simulink');
grid on;

figure;
plot(simout1.time, simout1.data + 2/3);
xlabel('Time (s)');
ylabel('Duty Ratio');
title('Change in FSFB duty ratio caused by step disturbance (Simulink)');
grid on;