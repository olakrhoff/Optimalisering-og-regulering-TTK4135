
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

C = [1 0];



[K,S,e] = dlqr(A,B,Q,R,[]);

eig(A-B*K);

p = (0.5) * [1; 1] + (0.03j) * [1; -1];

K_F = place(A', C', p);

phi = [A-B*K, B*K;
       diag([0,0]), A-K_F'*C];
eig(phi)
% Simulate
tf = 50; % Final time step
x = NaN(2,tf+1);
u = NaN(1,tf+1);
y = NaN(1,tf+1);
x_hat = NaN(2,tf+1);

x0     = [5, 1]'; % Initial state
x0_hat = [6, 0]'; % Initial state estimate

x(:,1) = x0;
x_hat(:,1) = x0_hat;

for t = 1:tf
    % System simulated one step ahead:
    u(:,t) = -K*x_hat(:,t);
    x(:,t+1) = A*x(:,t) + B*u(:,t);
    y(:,t) = C*x(:,t);
    % Calculate state estimate based on measurement y:
    x_hat(:,t+1) = A*x_hat(:,t) + B*u(:,t) + K_F'*(y(:,t) - C*x_hat(:,t)); %Slide 10: Lecture; MPC example: Adaptive Cruise Control
end

%% Plot
t_vec = 0:tf; % Time vector

% Plot optimal trajectory
figure(1);
subplot(2,1,1);
plot(t_vec, x, '--', 'linewidth', 2); hold on;
plot(t_vec, x_hat, '-'); hold off;
hleg = legend('$x_1(t)$', '$x_2(t)$', '$\hat{x}_1(t)$', '$\hat{x}_2(t)$');
set(hleg, 'Interpreter', 'Latex');
grid('on');
box('on');
ylim([-4, 8]);
ylabel('States and estimate');
subplot(2,1,2);
plot(t_vec,u);
box('on');
grid('on');
ylim([-8, 2]);
ylabel('u_t');
xlabel('t');