%Init known variables 
A = [1 0.5;
     0 1];
b = [0.125;
     0.5];
Q = [2 0;
     0 2];
R = 2;
N = 0;
 
%Solve with dlqr for P
[K,P,e] = dlqr(A,b,Q/2,R/2,N);

%Solve for the gain value K
R_inv = inv(R/2);
I = eye(2);
K = R_inv * b' * P * inv(I + b * R_inv * b' * P) * A;

%Check stability
eigen_values = eig(A-b*K);

%Check the largest value, if it holds, all other values hold
largest_eigen_value = max(abs(eigen_values));

