sys_tf = tf([3 2 1],[4 3 2 1]);

sys_ss = ss(sys_tf)

[A, B, C, D] = tf2ss(sys_tf.num{1}, sys_tf.den{1});

%cell array example
a = [1 2];
b = {(1), (2)};
b{2}(1);
c = {{1}, {'hello'}};
c{2}{1};

%bode plots are identical
%bode(sys_tf)
%hold
bode(sys_ss)