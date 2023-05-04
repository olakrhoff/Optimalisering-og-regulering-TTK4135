x_0 = [1.2, 1.2]';
x_0 = [-1.2, 1]';


[x_opt, x_iter, iter, alpha_rec] = steepest_descent(x_0);
%plot_iter_rosenbrock(x_iter);
[x_opt, x_iter, iter, alpha_rec] = newton(x_0);
%plot_iter_rosenbrock(x_iter);
[x_opt, x_iter, iter, alpha_rec] = bfgs(x_0);
%plot_iter_rosenbrock(x_iter);
x_opt
iter

%plot(alpha_rec)
%ylim([0 1.2])
%xlabel('Iteration') 
%ylabel('\alpha')
%title('\alpha-values at each iteration.')