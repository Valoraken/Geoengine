% Kinematic Measurement System Lab3
% Epoch 1
% Xiao Tianqi 3371477
% 
clear
close all
%% Data import and pre-processing
Data_1=xlsread('Messdaten_Trimble5601160512101415.687.xlsx');
Data_2=xlsread('Messdaten_Trimble5601160512101914.703.xlsx');

% Remove first and last second of each drive
% epoch_1
% using measurements 15- 315 327-600

T_1 = Data_1([52:315 327:547] , 1);             % time
Hz_1 = Data_1([52:315 327:547], 2);            % horizontal angle
V_1 = Data_1([52:315 327:547], 3);             % vertical angle
S_1 = Data_1([52:315 327:547], 4);             % slope distance
Y_1 = Data_1([52:315 327:547], 5);             % Y coordinates
X_1 = Data_1([52:315 327:547], 6);             % X coordinates
ds_1 = Data_1([52:315 327:547], 7);            % lateral deviation
n_1 = length(X_1);              % number of measurements
D_1=zeros(n_1, 1);              % distance
dt_1 = zeros(n_1, 1);           % time diiferences

% epoch_2
% using measurements 86-222 236-388
T_2 = Data_2([86:222 236:388], 1);
Hz_2 = Data_2([86:222 236:388], 2);
V_2 = Data_2([86:222 236:388], 3);
S_2 = Data_2([86:222 236:388], 4);
Y_2 = Data_2([86:222 236:388], 5);
X_2 = Data_2([86:222 236:388], 6);
ds_2 = Data_2([86:222 236:388], 7);
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
ds_mov1=[movmean(ds_1(1:301),5,'Endpoints','discard');movmean(ds_1(302:575),5,'Endpoints','discard')];
ds_mov2=[movmean(ds_2(1:137),5,'Endpoints','discard');movmean(ds_2(138:290),5,'Endpoints','discard')];
delta_ds1=max(ds_mov1)-min(ds_mov1);
delta_ds2=max(ds_mov2)-min(ds_mov2);
% Calculate velocity
v_1=[0;D_1(2:301)./(dt_1(2:301)/1000);D_1(303:575)./(dt_1(303:575)/1000)];
v_2=[0;D_2(2:137)./(dt_2(2:137)/1000);D_2(139:290)./(dt_2(139:290)/1000)];
vmean_1=sum(v_1)/length(v_1);
vmean_2=sum(v_2)/length(v_2);

%
delta_t1=delta_ds1/vmean_1;
delta_t2=delta_ds2/vmean_2;

% visualize
figure(1)
hold on
plot(X_1,Y_1);
figure(2)
plot(X_2,Y_2);
figure(3)
plot(1:n_1,X_2)