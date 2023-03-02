N = 30;
r = 1;

A = [0 0 0;
        0 0 1;
        0.1 -0.79 1.78];
B = [1; 0; 0.1];
C = [0 0 1];
R = 2*r;
Q = 2*diag([0 0 1]);

x_0 = [0; 0; 1];

% Now we have all the values needed to start solving the EQP

Q_n = kron(eye(N), Q);
R_n = kron(eye(N), R);
G = blkdiag(Q_n, R_n);

x_count = size(A,2); % Columns in A
u_count = size(B,2); % Columns in B

c = zeros(N*(x_count + u_count), 1); % We do not have a linear part to our model, hence just zeros
zero = zeros(N*x_count);

Aeq_c1 = sparse(eye(N*x_count));                         % Component 1 of A_eq
Aeq_c2 = sparse(kron(diag(ones(N-1,1),-1), -A));   % Component 2 of A_eq
Aeq_c3 = sparse(kron(eye(N), -B));                    % Component 3 of A_eq
A_eq = [Aeq_c1 + Aeq_c2, Aeq_c3];

%A_eq = [eye(N*x_count) + kron(diag(ones(N-1, -1), -A)), kron(eye(N), -B)]
b_eq = [A*x_0; zeros((N-1)*x_count, 1)];

solve = [G -A_eq'; A_eq zero]\[-c; b_eq];

% Extracting variables
z = solve(1:N*(x_count+u_count));   % Variable vector (solve includes lambdas)
y = [x0(3); z(x_count:x_count:N*x_count)]; % y = x3
u = z(N*x_count+1:N*x_count+N*u_count);    % Control
% Time vector
t = 1:N;

% Plot optimal trajectory
figure(1);
subplot(2,1,1);
plot([0,t],y,'-ko'); % Plot on 0 to N
grid('on');
ylabel('y_t')
subplot(2,1,2);
plot(t-1,u,'-ko'); % Plot on 0 to N-1
grid('on');
xlabel('t');
ylabel('u_t');
