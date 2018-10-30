d1 = [1 5 4 0];
L1 = tf(4,d1);
T1 = L1/(1+L1);
step(T1)
hold on
L2 = tf(5,d1);
T2 = L2/(1+L2);
step(T2)
L3 = tf(16,d1);
T3 = L3/(1+L3);
step(T3)
legend('K = 4','K = 5','K = 16')