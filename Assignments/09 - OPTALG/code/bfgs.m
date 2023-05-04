function [x_opt, x_iter, iter, alpha_rec] = bfgs(x_0)

max_itr = 2000;
grad_tol = 0.001;

n = size(x_0,1);

alpha_rec = NaN(1,max_itr);
x_iter = NaN(n,max_itr);

k = 1;
x = x_0;
x_iter(:, k) = x;
fval = f(x);
grad = gradient(x);
H = eye(n);

while (k < max_itr) && (norm(grad) >= grad_tol)
    fval = f(x);
    p = -H*grad;
    alpha_0 = 1;
    alpha = linesearch(x, p, fval, grad, alpha_0);
    alpha_rec(k) = alpha;
    x_prev = x;
    x = x + alpha*p;
    x_iter(:, k+1) = x;
    grad_prev = grad;
    grad = gradient(x);
    
    s_k = x - x_prev;
    y_k = grad - grad_prev;
    H = H_kp1(H, s_k, y_k);
    k = k+1;
end
fval = f(x);

% Return values
x_opt = x;
fval_opt = f(x_opt);
f_iter = fval;
iter = k;
end

function p = newton_func(x, grad)
    B = hessian(x);
    p = -B\grad; % inv(B)*grad
end

function alpha_k = linesearch(xk, pk, fk, gradk, alpha_0)
    alpha = alpha_0;
    rho = 0.95;
    c1 = 1e-4;
    while f(xk + alpha*pk) > fk + c1*alpha*gradk'*pk
        alpha = rho*alpha;
    end
    alpha_k = alpha;
end

function fval = f(x)
    fval = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
end

function grad = gradient(x)
    grad = [ -400*(x(1)*x(2)-x(1)^3) + 2*x(1) - 2 ;
             200*(x(2)-x(1)^2)                   ];
end

function H_next = H_kp1(H_k, s_k, y_k)
    rho_k = 1/(y_k'*s_k);
    I = eye(numel(s_k));
    H_next = (I - rho_k*s_k*y_k')*H_k*(I - rho_k*y_k*s_k') + rho_k*(s_k*s_k');
end