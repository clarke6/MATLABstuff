n1 = [1 1];
d1 = [1 0 0];
T1 = tf(n1, d1);
rlocus(T1)
figure
n2 = 1;
d2 = [1 4 4 0];
T2 = tf(n2, d2);
rlocus(T2)
figure
d3 = [1 20 101 0];
T3 = tf(n2, d3);
rlocus(T3)  