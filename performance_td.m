function z = perfomance_td(x,des_perfomance,weights)
% x(1) = ki , x(2) = ka

% Optimal Solution given in the paper
% ka = -0.006;
% ki = 8.61;

kr = 0.2903 ; 
s = tf('s');
wd = 22.4;
ksi_d = 0.052;
K1 = -1116.5;
K3 = 0.6477;
T_alpha = 0.676;
A11 = 0.001054;
A12 = -0.00081;

G_q = K3*(1+T_alpha*s)*wd^2/(s^2+2*ksi_d*wd*s+wd^2);

G_z = K1*(1+A11*s+A12*s^2)*wd^2/(s^2+2*ksi_d*wd*s + wd^2);

ksi_a = 0.7;
wa = 250;
ksi_r = 0.65;
wr = 500;

G_acc = wa^2/(s^2+2*ksi_a*wa*s+wa^2);

G_gyro = wr^2/(s^2+2*ksi_r*wr*s+wr^2);



w1 = weights.w1; % weight for undershoot error
w2 = weights.w2; % weight for overshoot error
w3 = weights.w3; % weight for settling time error
w4 = weights.w4; % weight for risetime error

usd = des_perfomance.us; % desired us
osd = des_perfomance.os; % desired os
tsd = des_perfomance.ts; % desired ts
rsd = des_perfomance.rs; % desired rs

G_ol = (x(2)*x(1)*(1/s)*kr*G_acc*G_z)/(1 + kr*G_acc*G_gyro*G_q + kr*G_acc*G_q*G_gyro*x(1)*(1/s));


% G_ol = minreal(G_ol);

G_cl = G_ol/(1+G_ol);

G_cl = minreal(G_cl); 

y = stepinfo(G_cl); % gives the time domain perdomance criterion

us = y.Undershoot; % System undershoot
os = y.Overshoot;  % System overshoot
ts = y.SettlingTime; % System settlingtime
rs = y.RiseTime; % System risetime

z = w1*(us-usd)^2 + w2*(os-osd)^2 + w3*(ts-tsd)^2 + w4*(rs-rsd)^2;


end