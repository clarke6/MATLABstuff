clear
close all
format compact

a = [0,1;0,0]
b = [0;1]
c = [1,0]
d = 0

sys_ss = ss(a,b,c,d)
sys_tf = tf(sys_ss)

p = [-2, -5]
k = acker(a,b,p)
a_cl = a - b*k
b_cl = b
c_cl = c
d_cl = d

figure(1)
step(ss(a_cl,b_cl,c_cl,d_cl))

p = [-2, -5]
k = acker(a,b,p)
a_cl = a - b*k
b_cl = b
c_cl = c
d_cl = d

figure(2)
step(ss(a_cl,b_cl,c_cl,d_cl))

cx = eye(size(a));
dx = zeros(size(b));

Nbar = -1/(c*inv(a-b*k)*b)

figure(3)
step(ss(a_cl,b_cl*Nbar,c_cl,d_cl))

aa = [a [0;0];-c,0];
bb = [b;0];
cc = [-c,0];
pp = [-1+j, -1-j, -5];
ka = acker(aa,bb,pp)
