load('Time.mat')
load('observed_10_3.mat')
x1m = observed_10_3(2,:);              % State x1 from solution
x2m = observed_10_3(3,:);              % State x2 from solution
x3m = observed_10_3(4,:);              % State x3 from solution
x4m = observed_10_3(5,:);              % State x4 from solution

% figure();
% plot(t,x1,'r', Time,x1m,'b', 'LineWidth', 2); grid on;
% legend('\lambda^*','\lambda')
% 
% figure(); 
% plot(t, x2,'o-r',Time, x2m,'b', 'LineWidth', 2); grid on;
% legend('r^*','r')
% 
% figure();
% plot(t,x3,'o-r',Time, x3m, 'b','LineWidth', 2); grid on;
% legend('p^*','p')
% 
% figure();
% plot(t,x4,'o-r',Time, x4m,'b','LineWidth', 2); grid on;
% legend('p_{dot}^*','p_{dot}')
% 
% figure();
% plot(t,x5,'o-r',Time, x5m, 'b','LineWidth', 2); grid on;
% legend('e^*','e')
% 
% figure();
% plot(t,x6,'-r',Time, x6m,'b', 'LineWidth', 2); grid on;
% legend('e_{dot}^*','e_{dot}');


% Input imported to helicopter
x_opt = [x1 x2 x3 x4];
calculated_input.time = t;
calculated_input.signals.values = u;
calculated_input.signals.dimensions = 2;
figure();
plot(t, calculated_input.signals.values, '-o', 'Linewidth', 2); grid on;
legend({'p_c', 'e_c'}, 'location', 'NorthEast', 'FontSize', 36)

%% [---- Task 10_4_4 - With feedack ----]
close all

x1m_lqr = observed_10_3(2,:);              % State x1 from solution
x2m_lqr = observed_10_3(3,:);              % State x2 from solution
x3m_lqr = observed_10_3(4,:);              % State x3 from solution
x4m_lqr = observed_10_3(5,:);              % State x4 from solution

figure();
plot(t,x1,'r', Time(1,2400:end),x1m_lqr(2400:end),'b',Time,x1m,'--g', 'LineWidth', 2); grid on;
legend('\lambda^*','\lambda')
axis([5 20 -0.5 3.5])

figure(); 
plot(t(20:end), x2(20:end),'r',Time(1,2400:end), x2m_lqr(2400:end),'b', 'LineWidth', 2); grid on;
legend('r^*','r')
axis([5 20 -1 0.2])

figure();
plot(t(20:end),x3(20:end),'r',Time(1,2400:end), x3m_lqr(2400:end), 'b','LineWidth', 2); grid on;
legend('p^*','p')
axis([5 20 -0.6 0.6])

figure();
plot(t(20:end),x4(20:end),'r',Time(1,2400:end), x4m_lqr(2400:end),'b','LineWidth', 2); grid on;
legend('p_{dot}^*','p_{dot}')
axis([5 20 -1 0.8])