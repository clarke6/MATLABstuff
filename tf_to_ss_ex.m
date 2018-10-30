sys_tf=tf([2,2],[1,2]);
sys_ss=ss(-2,1,-2,2); %See notes for converting ts to ss
bode(sys_ss)
figure()
bode(sys_tf)