x0 = [1.2, 1.2]'; % initial point
x0 = [-1.2, 1]'; % initial point
[x, fval, x_iter, iter] =  nelder_mead(x0, 'report');
plot_iter_rosenbrock(x_iter);
iter