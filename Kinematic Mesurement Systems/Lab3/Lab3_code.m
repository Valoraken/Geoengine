% Kinematic Measurement System Lab3
% Xiao Tianqi 3371477
% 28.05.2019

clear
close all
clc
set(0,'defaultfigurecolor','w')
%% Data import and pre-processing
Data_1=xlsread('Messdaten_Trimble5601160512101415.687.xlsx');
Data_2=xlsread('Messdaten_Trimble5601160512101914.703.xlsx');
Xs=0;
Ys=0;
Ori=0;
% Remove first and last second of each drive
% epoch_1
% using measurements 10-600
T_1 = Data_1(10:600 , 1)/1000;       % time [s]
Hz_1 = Data_1(10:600, 2);            % horizontal angle [gon]
V_1 = Data_1(10:600, 3);             % vertical angle [gon]
S_1 = Data_1(10:600, 4);             % slope distance [m]
Y_1 = Data_1(10:600, 5);             % Y coordinates [m]
X_1 = Data_1(10:600, 6);             % X coordinates [m]
ds_1 = Data_1(10:600, 7);            % lateral deviation [m]
n_1 = length(X_1);              % number of measurements
D_1=zeros(n_1, 1);              % distance
dt_1 = zeros(n_1, 1);           % time diiferences

% epoch_2
% using measurements 86-388
T_2 = Data_2(86:388, 1)/1000;
Hz_2 = Data_2(86:388, 2);
V_2 = Data_2(86:388, 3);
S_2 = Data_2(86:388, 4);
Y_2 = Data_2(86:388, 5);
X_2 = Data_2(86:388, 6);
ds_2 = Data_2(86:388, 7);
n_2 = length(X_2);
D_2=zeros(n_2, 1);
dt_2 = zeros(n_2, 1);

% calculate distance and time differences
for i = 1:n_1-1
    D_1(i+1) = sqrt((X_1(i+1)-X_1(i))^2+(Y_1(i+1)-Y_1(i))^2);
    dt_1(i+1) = T_1(i+1)-T_1(i);
end

for i = 1:n_2-1
    D_2(i+1) = sqrt((X_2(i+1)-X_2(i))^2+(Y_2(i+1)-Y_2(i))^2);
    dt_2(i+1) = T_2(i+1)-T_2(i);
end
%% Calculation
% moving average of lateral deviation
ds_mov1=movmean(ds_1,5);
ds_mov2=movmean(ds_2,5);
ds1_mean=mean(ds_1);
ds2_mean=mean(ds_2);
delta_ds1=max(ds_mov1)-min(ds_mov1);
delta_ds2=max(ds_mov2)-min(ds_mov2);
% Calculate velocity
v_1=[0;D_1(2:length(D_1))./dt_1(2:length(dt_1))];
v_2=[0;D_2(2:length(D_2))./dt_2(2:length(dt_2))];
% mean velocity
vmean_1=mean(v_1);
vmean_2=mean(v_2);

% determine sychronization error
delta_t1=delta_ds1/vmean_1;
delta_t2=delta_ds2/vmean_2;

%% Coordinates Correction
% epoch 1
SH_1=S_1.*sin(V_1./200.*pi);
Sreal_1=zeros(size(SH_1));
for i=1:n_1-1
    Sreal_1(i+1)=SH_1(i+1)+delta_t1./dt_1(i+1).*(SH_1(i)-SH_1(i+1));
end

X1=Xs+Sreal_1.*cos(Ori+Hz_1./200.*pi);
Y1=Ys+Sreal_1.*sin(Ori+Hz_1./200.*pi);

% epoch 2
SH_2=S_2.*sin(V_2./200.*pi);
Sreal_2=zeros(size(SH_2));
for i=1:n_2-1
    Sreal_2(i+1)=SH_2(i+1)+delta_t1./dt_2(i+1).*(SH_2(i)-SH_2(i+1));
end

X2=Xs+Sreal_2.*cos(Ori+Hz_2./200.*pi);
Y2=Ys+Sreal_2.*sin(Ori+Hz_2./200.*pi);

%% visualize
% Coordinates of the prism
figure(1)
hold on
plot(X_1,Y_1);
xlabel('X[m]');
ylabel('Y[m]');
title('Coordinates of the prism in velocity 1')

figure(2)
hold on
plot(X_2,Y_2);
xlabel('X[m]');
ylabel('Y[m]');
title('Coordinates of the prism in velocity 2')

% Lateral Deviation and Smoothed Lateral Deviation
figure(3)
% plot(1:n_1,X_1)
hold on
plot(T_1-T_1(1),ds_1,'Color', [0.5 0.4 0.9]);
plot(T_1-T_1(1),ds_mov1,'r -- ');
xlabel('T [s]');
ylabel('d [m]');
legend('Lateral Deviation','Smoothed Lateral Deviation');
title('Lateral Deviation and Smoothed Lateral Deviation of velocity 1')

figure(4)
hold on
plot(T_2-T_2(1),ds_2,'Color', [0.9 0.7 0.4]);
plot(T_2-T_2(1),ds_mov2,'b -- ');
xlabel('T [s]');
ylabel('d [m]');
legend('Lateral Deviation','Smoothed Lateral Deviation');
title('Lateral Deviation and Smoothed Lateral Deviation of velocity 2')

% Velocity 
figure (5)
hold on
plot(T_1(2:size(T_1))-T_1(1),v_1(2:size(v_1)),'Color', [0.7 0.6 0.9]);
xlabel('T [s]');
ylabel('V [m/s]');
title('Velocity of Measurement 1')

figure (6)
hold on
plot(T_2(2:size(T_2))-T_2(1),v_2(2:size(v_2)),'Color', [0.9 0.2 0.6]);
xlabel('T [s]');
ylabel('V [m/s]');
title('Velocity of Measurement 2')

% Corrected Coordinates and Original Coordinates
figure (7)
hold on
plot(X_1,Y_1);
plot(X1(2:size(X1)),Y1(2:size(Y1)));
legend('Original Coordinates ','Corrected Coordinates');
xlabel('X[m]');
ylabel('Y[m]');
title('Corrected Coordinates and Original Coordinates of the prism in velocity 1')

figure (8)
hold on
plot(X_2,Y_2);
plot(X2(2:size(X2)),Y2(2:size(Y2)));
legend('Original Coordinates ','Corrected Coordinates');
xlabel('X[m]');
ylabel('Y[m]');
title('Corrected Coordinates and Original Coordinates of the prism in velocity 2')