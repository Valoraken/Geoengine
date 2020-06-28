% Lab4 Task 1
close all
clear
clc
[G, S, R, E, C, numEpoche, numSat,numSat_epoch] = readObs('INSA00DEU_R_20193391200_01H_01S_MO.rnx')
% load('Rinex.mat')
%% a) Interpretation of the header 
%   type of receiver and software version 
%       Type:		TRIMBLE NETR9       
%       Version:	5.37    
%   GNSS systems from which observations are included 
%       G : GPS
%       S: SBAS payload
%       R: Glona2ss
%       E: Gallieo
%       C: Beidou
%   starting epoch in the corresponding time system
%   	2019    12     5    12     0    0.0000000     GPS time        
%   data rate of the observations 
%     	1.000 observation/s
%   main types of observations and their units
%       C = Code / Pseudorange      [meters] 
%       L = Phase                   [full cycles]
%       S = Raw signal strength 	[dBHz]            (signal to noise ratio) 
%% b) Plot
figure()
hold on
t=1:size(numSat,2);
plot(t,numSat)
title(' total number of received satellites from all systems ')

figure()
hold on
plot(t,numSat_epoch)
title(' the number of received satellites separated by each system ')
legend('G','S','R','E','J','C');
%% c)
% signal to noise ratio -> plot S
% GPS PRN06
S_G=[G{1,6}.S1C';G{1,6}.S2W';G{1,6}.S2X';G{1,6}.S5X'];
% GLONASS PRN22
S_R=[R{1,22}.S1C';R{1,22}.S1P';R{1,22}.S2C';R{1,22}.S2P'];
figure()
plot(t,[S_G;S_R])
legend('GPS-S1C','GPS-S2W','GPS-S2X','GPS-S5X','GLONASS-S1C','GLONASS-S1P','GLONASS-S2C','GLONASS-S2P');

figure()
subplot(1,2,1)
hold on
plot(t,S_G)
legend('GPS-S1C','GPS-S2W','GPS-S2X','GPS-S5X')

subplot(1,2,2)
hold on
plot(t,S_R)
legend('GLONASS-S1C','GLONASS-S1P','GLONASS-S2C','GLONASS-S2P')

%% d)
% phase observations -> plot L

% GPS PRN06
L_G=[G{1,6}.L1C';G{1,6}.L2W';G{1,6}.L2X';G{1,6}.L5X'];
% GLONASS PRN22
L_R=[R{1,22}.L1C';R{1,22}.L1P';R{1,22}.L2C';R{1,22}.L2P'];

figure()
plot(t,[L_G;L_R])
legend('GPS-L1C','GPS-L2W','GPS-L2X','GPS-L5X','GLONASS-L1C','GLONASS-L1P','GLONASS-L2C','GLONASS-L2P');

figure()
subplot(1,2,1)
hold on
plot(t,L_G)
legend('GPS-L1C','GPS-L2W','GPS-L2X','GPS-L5X')

subplot(1,2,2)
hold on
plot(t,L_R)
legend('GLONASS-L1C','GLONASS-L1P','GLONASS-L2C','GLONASS-L2P')
%% e)
% i==2 j==5
% phi:L1C
% r = C1C
% GPS frequency
L1=1575.42*1e6;   
L2=1227.6*1e6;
L5=1176.45*1e6;
% Glonass frequency
k=-3
G1=(1602+k*9/16)*1e6;
G2=(1246+k*7/16)*1e6;

% GPS: 1C 2W 2X 5X
g1=CMC(L1,L2,G{1,6}.L1C,G{1,6}.L2W,G{1,6}.C1C,G{1,6}.C2W)';
g2=CMC(L2,L5,G{1,6}.L2W,G{1,6}.L5X,G{1,6}.C2W,G{1,6}.C5X)';
g3=CMC(L2,L5,G{1,6}.L2X,G{1,6}.L5X,G{1,6}.C2X,G{1,6}.C5X)';
g4=CMC(L1,L2,G{1,6}.L1C,G{1,6}.L2X,G{1,6}.C1C,G{1,6}.C2X)';
g5=CMC(L1,L5,G{1,6}.L1C,G{1,6}.L5X,G{1,6}.C1C,G{1,6}.C5X)';
% Glonass
r1=CMC(G1,G2,R{1,22}.L1C,R{1,22}.L2C,R{1,22}.C1C,R{1,22}.C2C)';
r2=CMC(G1,G2,R{1,22}.L1C,R{1,22}.L2P,R{1,22}.C1C,R{1,22}.C2P)';
r3=CMC(G1,G2,R{1,22}.L1P,R{1,22}.L2C,R{1,22}.C1P,R{1,22}.C2C)';
r4=CMC(G1,G2,R{1,22}.L1P,R{1,22}.L2P,R{1,22}.C1P,R{1,22}.C2P)';

% Plot
figure()
plot(t,[g1;g2;g3;g4;g5])
legend('1C-2W','2W-5X','2X-5X','1C-2X','1C-5X')
title('code minus carrier of GPS PRN06')
figure()
plot(t,[r1;r2;r3;r4])
legend('1C-2C','1C-2P','1P-2C','1P-2P')
title('code minus carrier of GLONASS PRN22')

%% Task 2 
% a)
TEC=132*1e16;
dL=40.3/L1^2*TEC;
TEC_L2=dL*L2^2/40.3/1e16;
% b)
d1=0.3;
d2=0.25;
d3=sqrt((L1^2/(L1^2-L2^2))^2*d1^2+(-L2^2/(L1^2-L2^2))^2*d2^2);

d1=0.25;
d5=0.2;
d6=sqrt((L1^2/(L1^2-L5^2))^2*d1^2+(-L5^2/(L1^2-L5^2))^2*d5^2);