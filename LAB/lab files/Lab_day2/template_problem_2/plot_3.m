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
ylabel("Pitch (rad)");
xlabel("Time (s)");
title("Plot of p and p_{ref}, with q = 1");
xlim([0,17]);


data_12 = load("data_12.mat");
data_1_2 = load("data_1_2.mat");
data_0_12 = load("data_0_12.mat");

plt_matrix = [data_12.data12'; data_1_2.data'; data_0_12.data'];

figure(6)
plot(plt_matrix(1,:),plt_matrix(2,:),plt_matrix(1,:),plt_matrix(4,:),plt_matrix(1,:),plt_matrix(6,:)),grid
legend("q = 12", "q = 1.2", "q = 0.12");
ylabel("u (rad)")
xlabel("Time (s)")
title("Comparison of input (u) for different q values");
xlim([0,25]);


%% travel comprasison
data_12 = load("travel_12.mat");
data_1_2 = load("travel_1.2.mat");
data_0_12 = load("travel_0.12.mat");

plt_matrix = [data_0_12.data'; data_1_2.data'; data_12.data'];
figure(7)
plot(plt_matrix(1,:),plt_matrix(2,:),plt_matrix(1,:),plt_matrix(3,:),plt_matrix(1,:),plt_matrix(4,:)),grid
legend("q = 12", "q = 1.2", "q = 0.12");
ylabel("Travel (rad)")
xlabel("Time (s)")
title("Comparison of estimated {\lambda} for different q values");
xlim([5,25]);


fys_data_12 = load("fys_travel_12.mat");
fys_data_1_2 = load("fys_travel_1.2.mat");
fys_data_0_12 = load("fys_travel_0.12.mat");


fys_plt_matrix = [fys_data_12.ans(1,1:141);fys_data_12.ans(2,1:141)+pi;data_12.data'];
figure(8)
plot(fys_plt_matrix(1,:),fys_plt_matrix(2,:),fys_plt_matrix(1,:),fys_plt_matrix(3,:)),grid
legend("{\lambda}", "{\lambda}*");
ylabel("Travel (rad)")
xlabel("Time (s)")
title("Comparison of observed {\lambda} and estimated {\lambda}*");
xlim([0,35]);


