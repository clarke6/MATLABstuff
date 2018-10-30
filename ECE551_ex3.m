%ctrb(A,B) : compute controllability matrix
A = [-4,1;0,-2]
B = [1;0]
CC = ctrb(A,B)
rank(CC)
AA = [-2,1;1,-2]
[v,d] = eig(AA)
norm(v(:,2),2)
norm(v(:,1),2)
a = [1,0;0,-3]
b = [1;1]
acker(a,b,[-1,-1])