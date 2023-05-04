function [x_opt, x_iter, iter, alpha_rec] = steepest_descent(x_0)

max_itr = 2000;
grad_tol = 0.001;

n = size(x_0,1);

alpha_rec = NaN(1,max_itr);
x_iter = NaN(n,max_itr);
% First iteration
k = 1;
%x(:,k) = x_0;
x = x_0;
x_iter(:, k) = x;
%fval(k) = f(x(:,k));
fval = f(x);
%grad(:,k) = gradient(x(:,k)); 
grad = gradient(x);
%p(:,k) = sd(grad(:,k));
p = sd(grad);
alpha_0 = 1;
%alpha(k) = linesearch(x(:,k), p(:,k), fval(k), grad(:,k), alpha_0);
alpha = linesearch(x, p, fval, grad, alpha_0);
alpha_rec(k) = alpha;
%x(:,k+1) = x(:,k) + alpha(k)*p(:,k);
x = x + alpha * p;
%grad(:,k+1) = gradient(x(:,k+1));
grad_prev = grad;
grad = gradient(x);
k = k+1;
x_iter(:, k) = x;

while (k < max_itr) && (norm(grad) >= grad_tol)
    fval = f(x);
    p_prev = p;
    p = sd(grad);
    alpha_0 = alpha * (grad_prev'*p_prev)/(grad'*p); % p. 59 N&W
    alpha = linesearch(x, p, fval, grad, alpha_0);
    alpha_rec(k) = alpha;
    x = x + alpha*p;
    x_iter(:, k+1) = x;
    grad_prev = grad;
    grad = gradient(x); % Done here since needed in the while condition
    k = k+1;
end
fval = f(x);

% Return values
x_opt = x;
fval_opt = f(x_opt);
f_iter = fval;
iter = k;
end

function p = sd(grad)
    p = -grad;
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