clear 
close all
format compact
s = tf('s');
sys_tf = 1/(s+1);
h=1/4;

gc = 4/s;

figure(1)
loop=sys_tf*h*gc;
margin(loop)

gc = 5.5/s;

figure(2)
loop=sys_tf*h*gc;
margin(loop)

sys_cl = sys_tf * gc/(1+loop);
figure(3)
step(sys_cl)