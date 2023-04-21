clear all; clc;

data=load('Excersize_5.mat');

array = data.ans;

t = array(1,:);
p_ref = array(2,:);
p = array(3,:);
travel = array(4,:);

figure(1);
plot(t, p_ref, t, p),grid;
legend("p_{ref}","p");
ylabel("Radians");
xlabel("Time (s)");
title("Plot of p and p_{ref}, with q = 1");
xlim([0,17]);


data_12 = load("template_problem_2\data_12.mat");
data_1_2 = load("template_problem_2\data_1_2.mat");
data_0_12 = load("template_problem_2\data_0_12.mat");

plt_matrix = [data_12.data12'; data_1_2.data'; data_0_12.data'];

figure(6)
plot(plt_matrix(1,:),plt_matrix(2,:),plt_matrix(1,:),plt_matrix(4,:),plt_matrix(1,:),plt_matrix(6,:)),grid
legend("q = 12", "q = 1.2", "q = 0.12");
ylabel("Radians")
xlabel("Time (s)")
title("Comparison of input (u) for different q values");
xlim([0,25]);