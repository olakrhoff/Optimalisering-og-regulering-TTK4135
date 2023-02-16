x1_1 = -1.5; x1_h = 6.5;
x2_1 = -1.5; x2_h = 6.5;
res = 0.01;
[x1, x2] = meshgrid(x1_1:res:x2_h, x1_1:res:x2_h);


f = -(3-0.4*x1).*x1 - (2-0.2*x2).*x2;
levels = (-12:2:8)';
[C,h] = contour(x1, x2, f, levels, 'Color', .7*[1 1 1]);
set(h, 'ShowText', 'on', 'LabelSpacing', 300);



G = [0.8 0;
    0 0.4];
c = [-3;-2];
A = [2 1;
    1 3];
b = [8; 15];



