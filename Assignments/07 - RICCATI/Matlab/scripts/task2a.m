
k1 = 1;
k2 = 1;
k3 = 1;
T = 0.1;
A = [1 T;
     -k2*T 1-k1*T];
B = [0;
     k3*T];
Q = [4 0;
     0 4];
R = 1;
N = 0;




[K,S,e] = dlqr(A,B,Q/2,R/2,N);
K

eig(A-B*K)