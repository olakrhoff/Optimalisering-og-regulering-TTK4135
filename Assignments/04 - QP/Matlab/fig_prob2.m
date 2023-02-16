x1_l = -1.5; x1_h = 6.5;
x2_l = -1.5; x2_h = 6.5;
res = 0.1;

[x1, x2] = meshgrid(x1_l:res:x2_h, x1_l:res:x2_h);

%f = -3*x1 - 2*x2;
f = -(3-0.4*x1).*x1 - (2-0.2*x2).*x2;
%levels = (-32:2:6)';
levels = (-12:2:8)';


fig = figure(1);
% Set fontsize for the labels and lengend
fontsize = 8;
% Set the figure dimension units to centimeters
set(gcf, 'units','centimeters');
% Figure paper dimensions:
w = 11; % 9 cm wide
h = 11; % 9 cm tall
% Set the dimensions of the figure
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1), pos(2), w, h]);

[C, h] = contour(x1, x2, f, levels, 'Color', .7*[1 1 1]);
% set(h, 'ShowText','on', 'LabelSpacing',500);
hold('on');
grid('off');
axis('square');

% plot([0 0], [x1_l, x1_h], 'k');
% plot([x2_l, x2_h], [0 0], 'k');

% plot_arrow(0,x1_l,0,x1_h);
point1_arrow = [0,x1_l];
point2_arrow1 = [0,x1_h];
d_point = point2_arrow1-point1_arrow;
quiver(point1_arrow(1), point1_arrow(2), d_point(1), d_point(2), 0, 'black')

% plot_arrow(x2_l,0,x2_h,0);
point1_arrow = [x2_l,0];
point2_arrow1 = [x2_h,0];
d_point = point2_arrow1-point1_arrow;
quiver(point1_arrow(1), point1_arrow(2), d_point(1), d_point(2), 0, 'black')


% Plot c_1 and c_2
x_1 = [-10,10];
x_2 = 8 - 2*x_1;
plot(x_1, x_2, 'k');
x_2 = (1/3)*(15 - x_1);
plot(x_1, x_2, 'k');

x1_opt = 1.8;
x2_opt = 4.4;
%plot(x1_opt, x2_opt, 'ok');
it1_x1 = 0;
it1_x2 = 0;
%plot(it1_x1, it1_x2, 'ok');
it2_x1 = 4;
it2_x2 = 0;
%plot(it2_x1, it2_x2, 'ok');

c1_x1 = -2;
c1_x2 = -1;

c2_x1 = -1;
c2_x2 = -3;

f_x1 = -3;
f_x2 = -2;

% plot_arrow(x1_opt,x2_opt, x1_opt+c1_x1,x2_opt+c1_x2, 'linewidth', 2);
point1_arrow = [x1_opt,x2_opt];
point2_arrow1 = [x1_opt+c1_x1,x2_opt+c1_x2];
d_point = point2_arrow1-point1_arrow;
%quiver(point1_arrow(1), point1_arrow(2), d_point(1), d_point(2), 0, 'black')

% plot_arrow(x1_opt,x2_opt, x1_opt+c2_x1,x2_opt+c2_x2, 'linewidth', 2);
point1_arrow = [x1_opt,x2_opt];
point2_arrow1 = [ x1_opt+c2_x1,x2_opt+c2_x2];
d_point = point2_arrow1-point1_arrow;
%quiver(point1_arrow(1), point1_arrow(2), d_point(1), d_point(2), 0, 'black')

% plot_arrow(x1_opt,x2_opt, x1_opt+f_x1,x2_opt+f_x2, 'linewidth', 2);
point1_arrow = [x1_opt,x2_opt];
point2_arrow1 = [ x1_opt+f_x1,x2_opt+f_x2];
d_point = point2_arrow1-point1_arrow;
%quiver(point1_arrow(1), point1_arrow(2), d_point(1), d_point(2), 0, 'black')


text_offset = 0.1;
%text(x1_opt+c1_x1-5*text_offset,x2_opt+c1_x2+0*text_offset, '$\nabla c_1$', 'fontsize', fontsize,'interpreter','latex')
%text(x1_opt+c2_x1-3*text_offset,x2_opt+c2_x2-2*text_offset, '$\nabla c_2$', 'fontsize', fontsize,'interpreter','latex')
%text(x1_opt+f_x1-1*text_offset,x2_opt+f_x2-2*text_offset, '$\nabla f$', 'fontsize', fontsize,'interpreter','latex')

% Iterations
%text(it1_x1+1*text_offset,it1_x2+2*text_offset, 'Iteration 1', 'fontsize', fontsize)
%text(it2_x1+1*text_offset,it2_x2+2*text_offset, 'Iteration 2', 'fontsize', fontsize)
%text(x1_opt+1*text_offset,x2_opt+2*text_offset, 'Iteration 3', 'fontsize', fontsize)

xlabel('x_1');
ylabel('x_2');
% Set the fontsize of all elements
set(gca, 'fontsize', fontsize); % Sets fontsize on legend and numbers on axis
set(get(gca, 'XLabel'), 'fontsize', fontsize); % Sets the fontsize for the xlabel
set(get(gca, 'YLabel'), 'fontsize', fontsize); % Sets the fontsize for the ylabel

set(gca,'XTick',[-1:1:6],'YTick',[-1:1:6])
xlim([-1.5,6.5])
ylim([-1.5,6.5])

hold('off');

% Create the pdf file
% pdfmatlabfrag(gcf, 'fig_prob2.pdf');