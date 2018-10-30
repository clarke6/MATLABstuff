close all


L = 560e-6;
C = 100e-6;
R = 5;

sys_ss = ss([-1/(R*C),1/C;-1/L,0],[0;1/L],[1,0],0);
sys_tf = tf([1/(L*C)],[1, 1/(R*C), 1/(L*C)]);
sys_ssPV = ss([0,1;-1/(L*C),-1/(R*C)],[0;1],[1/(L*C),0],0);
sys_ssDPV = ss([-1/(R*C),1;-1/(L*C),0],[0;1/(L*C)],[1,0],0);
bode(sys_tf)
title('Transfer Function')
figure()
bode(sys_ssPV)
title('Phase Variable')
figure()
bode(sys_ssDPV)
title('Dual Phase Variable')
figure()
bode(sys_ss)
title('Original State Space Description')
figure()
bode(sys_ss)
hold on
bode(sys_tf)
bode(sys_ssPV)
bode(sys_ssDPV)
title('All Plots overlapping')