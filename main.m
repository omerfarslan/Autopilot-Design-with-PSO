% Omer Faruk ARSLAN
% 2231256

close all;
clear;
clc;

rad = pi/180;

kr = 0.2903 ; % chosen by looking bandwidth ratio(1/3)
A =  10^(5/20);   %1.9953; %  6-dB gain margin(required, design parameter)
theta = 30*rad; % phase margin must be radian


% 30 deg & 5 dB works fine

%% Vehicle Dynamics & TFs

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

%% Actuator & Rate Gyro TFs

ksi_a = 0.7;
wa = 250;
ksi_r = 0.65;
wr = 500;

G_acc = wa^2/(s^2+2*ksi_a*wa*s+wa^2);

G_gyro = wr^2/(s^2+2*ksi_r*wr*s+wr^2);


%% FIRST PHASE => STABILITY

logg = linspace(1,100,200);


xsoln1 =[]; % inner
xsoln2 =[]; % inter
xsoln3 =[]; % outer

for i=1:length(logg)

    w = logg(i);

    % Modified Version to Handle Large Number effects on numerical solution
    Ni = 1;
    Di = j*w;
    Nacc = wa^2;
    Dacc = ((j*w)^2+2*ksi_a*wa*j*w+wa^2);
    Ngyro = wr^2;
    Dgyro = ((j*w)^2+2*ksi_r*wr*j*w+wr^2);
    Nq = K3*(1+T_alpha*j*w);
    Dq = ((j*w)^2+2*ksi_d*wd*j*w+wd^2)/wd^2;
    Nz = K1*(1+A11*(j*w)+A12*(j*w)^2);
    Dz = ((j*w)^2+2*ksi_d*wd*(j*w)+wd^2)/wd^2;


    % for outer loop

    P0 = kr*Nacc*Nq*Ngyro*Ni*Dz;
    Q0 = A*exp(-j*theta)*kr*Ni*Nacc*Nz*Dq*Dgyro;
    R0 = Di*Dacc*Dz*Dq*Dgyro + kr*Nacc*Ngyro*Nq*Di*Dz;

    P_real = double(real(P0)); % to solve with fsolve function
    P_imag = double(imag(P0));
    Q_real = double(real(Q0));
    Q_imag = double(imag(Q0));
    R_real = double(real(R0));
    R_imag = double(imag(R0));

    x0 = [1;-1]; %we know ka will be negative


    F = @(x) [P_real*x(1)+Q_real*x(1)*x(2)+R_real;
        P_imag*x(1)+Q_imag*x(1)*x(2)+R_imag];


    options = optimoptions('fsolve','MaxFunEvals',500,'Tolfun',1e-3);

    [x,fval,exitflag] = fsolve(F,x0,options);

    if exitflag > 0
        xsoln3 = [xsoln3 x];
    end

    clear x;

    % intermediate loop

    P0 = A*exp(-j*theta)*kr*Nacc*Nq*Ngyro*Ni*Dz;
    Q0 = A*exp(-j*theta)*kr*Nacc*Ni*Nz*Dq*Dgyro;
    R0 = Dacc*Di*Dq*Dgyro*Dz + kr*Nacc*Ngyro*Nq*Di*Dz;


    P_real = double(real(P0));
    P_imag = double(imag(P0));
    Q_real = double(real(Q0));
    Q_imag = double(imag(Q0));
    R_real = double(real(R0));
    R_imag = double(imag(R0));

    x0 = [1;-1]; %we know ka will be negative


    F = @(x) [P_real*x(1)+Q_real*x(1)*x(2)+R_real;
        P_imag*x(1)+Q_imag*x(1)*x(2)+R_imag];


    options = optimoptions('fsolve','MaxFunEvals',500,'Tolfun',1e-3);

    [x,fval,exitflag] = fsolve(F,x0,options);

    if exitflag > 0
        xsoln2 = [xsoln2 x];
    end

    clear x;

    %     inner loop

    P0 = A*exp(-j*theta)*kr*Nacc*Nq*Ngyro*Ni*Dz;
    Q0 = A*exp(-j*theta)*kr*Ni*Nacc*Nz*Dq*Dgyro;
    R0 = Dacc*Di*Dq*Dgyro*Dz + A*exp(-j*theta)*kr*Nacc*Ngyro*Nq*Di*Dz;

    P_real = double(real(P0));
    P_imag = double(imag(P0));
    Q_real = double(real(Q0));
    Q_imag = double(imag(Q0));
    R_real = double(real(R0));
    R_imag = double(imag(R0));

    x0 = [1;-1]; %we know ka will be negative


    F = @(x) [P_real*x(1)+Q_real*x(1)*x(2)+R_real;
        P_imag*x(1)+Q_imag*x(1)*x(2)+R_imag];

    options = optimoptions('fsolve','MaxFunEvals',500,'Tolfun',1e-3);


    [x,fval,exitflag] = fsolve(F,x0,options);


    if exitflag > 0
        xsoln1 = [xsoln1 x];

    end


    clear x;

end

xq1 = [-0.05:0.001:0];
vq1 = interp1(xsoln1(2,1:(end-18)),xsoln1(1,1:(end-18)),xq1,'linear','extrap');

xq2 = [-0.05:0.001:0];
vq2 = interp1(xsoln2(2,:),xsoln2(1,:),xq2,'spline','extrap');

xq3 = [-0.05:0.001:0];
vq3 = interp1(xsoln3(2,:),xsoln3(1,:),xq3,'linear','extrap');


figure;
plot(xq1,vq1,'LineWidth',2);
hold on;
plot(xq2,vq2,'LineWidth',2);
hold on;
plot(xq3,vq3,'LineWidth',2);
title('Gain Boundaries for Given Phase & Gain Margin');
xlabel('acc gain ka');
ylabel('integral gain ki');
legend('inner loop','inter loop','outer loop');
ylim([0 35]);
xlim([-0.05 0]);


%% SECOND PHASE => PSO SEARCH REGION

k_al = -0.02;
k_au = 0;
k_il = 0;

for i = 1:length(xq1)

    k_iu(i) = min([vq1(i) vq2(i) vq3(i)]);

end
ka_vec = [-0.02:0.0001:0];

% plot(xq1,k_iu);
figure;
area(xq1,k_iu);
title('Search Region of PSO');
xlabel('Acc Gain Ka');
ylabel('Integral Gain Ki');


p = polyfit(xq1,k_iu,9);
y = polyval(p,ka_vec);

figure;
plot(ka_vec,y,'LineWidth',2);
title('Fitted Polynomial For PSO');
xlabel('AccGain Ka');
ylabel('Integral Gain Ki');

p_fun = @(x) p(1)*x^9 + p(2)*x^8 + p(3)*x^7 + p(4)*x^6 + p(5)*x^5 ...
    + p(6)*x^4 + p(7)*x^3 + p(8)*x^2 + p(9)*x + p(10);



%% THIRD PHASE => PERFORMANCE OPTIMIZATION with PSO

% use the notation of xsoln = [ki ka];

% Performance Criterion

weights.w1 = 1; % weight for undershoot error
weights.w2 = 1; % weight for overshoot error
weights.w3 = 10; % weight for settling time error
weights.w4 = 10; % weight for risetime error

% Desired Perfomance Respectivley
des_perfomance.us = 10;
des_perfomance.os = 2;
des_perfomance.ts = 0.3;
des_perfomance.rs = 0.15;


% Problem Definition

problem.CostFunction = @(x,des_performance,weights) performance_td(x,des_performance,weights);% cost function definition

problem.nVar = 2; % number of variables

problem.VarMinf =  [0 -0.02] ;  % lower bound for variables

problem.VarMaxf = @(x) [p_fun x];  % upper bound for variables

% PSO Paramaters

param.MaxIt = 100;
param.nPop = 10;
param.w = 1;
param.wdamp = 0.99;
param.c1 = 2;
param.c2 = 2;

% Call PSO % Results

disp('#############');
disp('PSO START');

out = PSO(problem,param,des_perfomance,weights);

disp('PSO END');
disp('#############');

seq = (param.MaxIt/5);
cost1_seq = out.BestCost_Array(2:seq);
cost2_seq = out.BestCost_Array(seq+1:2*seq);
cost3_seq = out.BestCost_Array(2*seq+1:3*seq);
cost4_seq = out.BestCost_Array(3*seq+1:4*seq);
cost5_seq = out.BestCost_Array(4*seq+1:5*seq);

figure;
plot(out.BestCost_Array(2:end),'LineWidth',2); % since first element is too high(assigned by me)
hold on;
scatter(seq,cost2_seq(1),100,'red','filled');
hold on;
scatter(2*seq,cost3_seq(1),100,'green','filled');
hold on;
scatter(3*seq,cost4_seq(1),100,'blue','filled');
hold on;
scatter(4*seq,cost5_seq(1),100,'black','filled');
title('PSO Cost Calculation');
xlabel('Iteration Number');
ylabel('Cost Value');
grid on;
legend('Cost over Time','PSO 1st RESET','PSO 2nd RESET','PSO 3rd RESET','PSO 4th RESET');

xcost = out.BestCost;
xsoln = double(out.BestSoln.Position);

disp(' ');
disp(['Best Solution = ' num2str(xsoln)]);
disp(['Best Cost = ' num2str(xcost)]);
disp('');

%% FINAL RESULTS & PLOTS

kr = 0.2903 ; % chosen by looking bandwidth ratio(1/3)

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

G_ol = (xsoln(2)*xsoln(1)*(1/s)*kr*G_acc*G_z)/(1 + kr*G_acc*G_gyro*G_q + kr*G_acc*G_q*G_gyro*xsoln(1)*(1/s));

G_cl = minreal(G_ol/(1+G_ol));

y = stepinfo(G_cl); % we can obtain required time domain data from here

% Compute Percent Errors

error_os = abs(100*(des_perfomance.os - y.Overshoot)/des_perfomance.os);
error_us = abs(100*(des_perfomance.us - y.Undershoot)/des_perfomance.us);
error_rs = abs(100*(des_perfomance.rs - y.RiseTime)/des_perfomance.rs);
error_ts = abs(100*(des_perfomance.ts - y.SettlingTime)/des_perfomance.ts);


disp('Time Domain Performance');
disp(['Error in Desired Overshoot    = %' num2str(error_os)]);
disp(['Error in Desired Undershoot   = %' num2str(error_us)]);
disp(['Error in Desired RiseTime     = %' num2str(error_rs)]);
disp(['Error in Desired SettlingTime = %' num2str(error_ts)]);

% After Time Domain Performance optimization, we can check frequency domain
% stability requirements (chosen as 5 dB & 30 deg at least)

G_ol_outer = (xsoln(2)*xsoln(1)*(1/s)*kr*G_acc*G_z)/(1 + kr*G_acc*G_gyro*G_q + kr*G_acc*G_q*G_gyro*xsoln(1)*(1/s));

G_ol_inter = (xsoln(2)*xsoln(1)*(1/s)*kr*G_acc*G_z + kr*G_acc*G_q*G_gyro*xsoln(1)*(1/s))/(1+kr*G_acc*G_gyro*G_q);

G_ol_inner = xsoln(2)*xsoln(1)*(1/s)*kr*G_acc*G_z + kr*G_acc*G_q*G_gyro*xsoln(1)*(1/s) + kr*G_acc*G_gyro*G_q;

% We can see that all 3 loops are satisfying the "stability requirements".

figure;
margin(G_ol_outer);
figure;
margin(G_ol_inter);
figure;
margin(G_ol_inner);

clc; % Write results clearly into command window

disp('#############');
disp(['Best Solution = ' num2str(xsoln)]);
disp(['Best Cost = ' num2str(xcost)]);
disp('');
disp('Time Domain Performance');
disp(['Error in Desired Overshoot    = %' num2str(error_os)]);
disp(['Error in Desired Undershoot   = %' num2str(error_us)]);
disp(['Error in Desired RiseTime     = %' num2str(error_rs)]);
disp(['Error in Desired SettlingTime = %' num2str(error_ts)]);
disp('');
disp('All 3 Loop has at least 5 dB GM and 30 deg PM.');
disp('#############');
