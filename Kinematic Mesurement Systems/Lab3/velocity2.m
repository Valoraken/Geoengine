% Kinematic Measurement System Lab3
% Kinematic Tachymeter Measurements - velocity 2
% Name: YI HONG
% Matriculation number:3294211
%
clc
clear
Data_v2 = load('Messdaten_Trimble5601180608134205_625-velocity2.xls');
T = Data_v2(:, 1);
Hz = Data_v2(:, 2);
V = Data_v2(:, 3);
S = Data_v2(:, 4);
Y = Data_v2(:, 5);
X = Data_v2(:, 6);
ds = Data_v2(:, 7);
n = length(X);
d = zeros(n, 1);
dt = zeros(n, 1);
for i = 1:n-1
    d(i+1) = sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2);
    dt(i+1) = T(i+1)-T(i);
end
v = d./(dt/1000);

v1 = v(45:202,:);    % preserve velocity at the position: 45-202 and 215-356
v2 = v(215:356,:);   % 1 represents the forward drive;2 represents the backward drive 
[p1, v1] = Sigma3(v1);  % 3 sigma principle
[p2, v2] = Sigma3(v2);            
v = [v1;v2];

T1 = T(44:202,:);
T1(p1+1) = [];
T2 = T(214:356,:);
T2(p2+1) = [];

Hz1 = Hz(44:202,:);
Hz1(p1+1) = [];
Hz2 = Hz(214:356,:);
Hz2(p2+1) = [];

V1 = V(44:202,:);
V1(p1+1) = [];
V2 = V(214:356,:);
V2(p2+1) = [];

S1 = S(44:202,:);
S1(p1+1) = [];
S2 = S(214:356,:);
S2(p2+1) = [];

X1 = X(44:202,:);
X1(p1+1) = [];
X2 = X(214:356,:);
X2(p2+1) = [];
X = [X1;X2];

Y1 = Y(44:202,:);
Y1(p1+1) = [];
Y2 = Y(214:356,:);
Y2(p2+1) = [];
Y = [Y1;Y2];

ds1 = ds(44:202,:);
ds1(p1) = [];
ds2 = ds(214:356,:);
ds2(p2) = [];
ds = [ds1;ds2];

%% Computation Process
% Smoothing of lateral deviation
ds_smooth1 = movmean(ds1, 5);  
ds_smooth2 = movmean(ds2, 5);
ds_smooth = movmean(ds, 5);

T_smooth1 = movmean(T1, 5);                                                % for plot purpose
T_smooth2 = movmean(T2, 5);                                                % for plot purpose

% Difference between max and min of lateral deviation
max_ds1 = max(ds_smooth1);
min_ds1 = min(ds_smooth1);
delta_ds1 = max_ds1 - min_ds1;

max_ds2 = max(ds_smooth2);
min_ds2 = min(ds_smooth2);
delta_ds2 = max_ds2 - min_ds2;

max_ds = max(ds_smooth);
min_ds = min(ds_smooth);
delta_ds = max_ds - min_ds;

% Mean velocity
v_mean1_new = mean(v1);
v_mean2_new = mean(v2);
v_mean_new = mean(v);

% Synchronisation error
detat1 = delta_ds1/v_mean1_new;
detat2 = delta_ds2/v_mean2_new;
detat = delta_ds/v_mean_new;

% Corrected coordinates
[X_corr1, Y_corr1] = CorrectedCoord(T1, Hz1, V1, S1, detat1);
[X_corr2, Y_corr2] = CorrectedCoord(T2, Hz2, V2, S2, detat2);

% Figures
figure(1);
plot(X, Y);
title('Coordinates of the prism (second drive)');
xlabel('X [m]');
ylabel('Y [m]');

figure(2);                                                                 % lateral deviations
plot(T1/1000,ds1, 'color', [0 0.5 1]);   
hold on;                                 
plot(T_smooth1/1000, ds_smooth1, '--', 'color', [1 0.5 0.5],'Linewidth', 2);
hold on;
plot(T2/1000,ds2, 'color', [0 0.5 1]);   
hold on;                                 
plot(T_smooth2/1000, ds_smooth2, '--', 'color', [1 0.5 0.5],'Linewidth', 2);
title('Lateral deviations and smoothed lateral devitations');
legend('Lateral deviations', 'Smoothed lateral devitations');
xlabel('Time [s]');
ylabel('Lateral deviation [m]');

figure(3);                
T1(1) = [];                                                                % delete time without velocity
plot(T1/1000, v1,'color', [0 0.5 1]);
hold on;
plot([T1(1)/1000, T1(end)/1000],[v_mean1_new, v_mean1_new],'--', 'color', [1 0.5 0.5],'Linewidth', 2);
hold on;
T2(1) = [];
plot(T2/1000, v2,'color', [0 0.5 1]);
hold on;
plot([T2(1)/1000, T2(end)/1000],[v_mean2_new, v_mean2_new],'--', 'color', [1 0.5 0.5],'Linewidth', 2);
legend('Velocity', 'mean velocity');
title('Velocity varies with time')
xlabel('Time [s]');
ylabel('Velocity [m/s]');

figure(4);
plot(X, Y, 'color', [0 0.5 1]);
hold on;
plot(X_corr1, Y_corr1, '--', 'color', [1 0.5 0.5],'Linewidth', 2);
hold on;
plot(X_corr2, Y_corr2, '--', 'color', [1 0.5 0.5],'Linewidth', 2);
title('Coordinates of the prism')
legend('Measured coordinate', 'Corrected coordinate');
xlabel('X [m]');
ylabel('Y [m]');