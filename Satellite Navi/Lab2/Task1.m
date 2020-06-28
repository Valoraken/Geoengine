close all
clear

%% Keplerian Elements
a = 26371000;      % semi-major axis [m]
e = 0.4;           % eccentricity
I = 45;            % inclination [deg]
omega = 110;       % argument of perigee [deg]
Omega  = 50;       % right ascension of the ascending node [deg]
t0=0;              % starting ephoch
M0=0               % mean anormaly [deg]
dt=10;             % sampling rate [s]

%% Compute positions and velocities for 1 day
T = 1*24*3600;     % Time Interval [s] (one day)
t = t0:dt:T;       % Time of each epoch [s]

[Position,Velocity] = kep2cart(a,e,deg2rad(I),deg2rad(omega),deg2rad(Omega),M0,t0,t);

figure
view (3)
hold on
grid on
scatter3(0,0,0,250,'b','filled')
plot3(Position(1,:),Position(2,:),Position(3,:))
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Orbit in geocentric coordinate system')

figure (2)
plot(t,Velocity);
legend('Velosity-x','Velosity-y','Velosity-z')
xlabel('Time [s]');
ylabel('Velosity [m/s]');
title('Velocity components')

%% Compute positions and velocities into five Keplerian Elements (a, e, I, ?, ¦Ø).
[a_c,e_c,I_c,W_c,w_c] = cart2kep(Position , Velocity);
% da=a_c-a;
% de=e_c-e;
% dI=rad2deg(I_c)-I;
% dW=rad2deg(W_c)-Omega;
% dw=rad2deg(w_c)-omega;

%% Transform the satellite orbit into the ECEF coordinate system

omega_E=2*pi/T; % angular velocity of earth [rad/s]
rot=angle2dcm(omega_E*t, zeros(size(t)), zeros(size(t)));
re =cell2mat( arrayfun(@(p) rot(:,:,p)*Position(:,p), 1:length(t), 'UniformOutput', false));

figure
view (3)
hold on
grid on
scatter3(0,0,0,600,'b','filled')
plot3(re(1,:),re(2,:),re(3,:))
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('Orbit in ECEF coordinate system')
%% Load Satellite orbit informations
GPS=importdata('GPS.txt',' ',1); % need Column 2-7
Galileo=importdata('Galileo.txt'," ",1); % need Column 3-8
GLONASS=importdata('GLONASS.txt'," ",1); % need Column 2-7
%% Compute Coordinates
% Difine Epochs
delta_t=30;
t = t0:delta_t:T;
% GPS data
[Position_GPS,Velocity_GPS] = kep2cart(1000*GPS.data(:,1),GPS.data(:,2),deg2rad(GPS.data(:,3)),deg2rad(GPS.data(:,5)),deg2rad(GPS.data(:,4)),deg2rad(GPS.data(:,6)),t0,t);
rot_GPS=angle2dcm(omega_E*t, zeros(size(t)), zeros(size(t)));
% GLONASS data
[Position_GLONASS,Velocity_GLONASS] = kep2cart(1000*GLONASS.data(:,2),GLONASS.data(:,3),deg2rad(GLONASS.data(:,4)),deg2rad(GLONASS.data(:,6)),deg2rad(GLONASS.data(:,5)),deg2rad(GLONASS.data(:,7)),t0,t);
rot_GLONASS=angle2dcm(omega_E*t, zeros(size(t)), zeros(size(t)));
% Galileo data
[Position_Galileo,Velocity_Galileo] = kep2cart(1000*Galileo.data(:,3),Galileo.data(:,4),deg2rad(Galileo.data(:,5)),deg2rad(Galileo.data(:,7)),deg2rad(Galileo.data(:,6)),deg2rad(Galileo.data(:,8)),t0,t);
rot_Galileo=angle2dcm(omega_E*t, zeros(size(t)), zeros(size(t)));
%% Plot Kepler Orbits for 1 day (sampling rate ?t = 30 s) of the GNSS satellites:
figure
view (3)
hold on
grid on
for i=1:length(GPS.data)
    Orbit=Position_GPS(:,(i-1)*size(t,2)+1:i*size(t,2));
    plot3(Orbit(1,:),Orbit(2,:),Orbit(3,:));
end
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
scatter3(0,0,0,250,'b','filled');
title('Orbit of GPS in geocentric coordinate system')

% figure
% view (3)
% hold on
% grid on
% for i=1:length(GPS.data)
%     Orbit=Position_GPS(:,(i-1)*size(t,2)+1:i*size(t,2));
%     re_GPS =cell2mat( arrayfun(@(p) rot_GPS(:,:,p)*Orbit(:,p), 1:length(t), 'UniformOutput', false));
%     plot3(re_GPS(1,:),re_GPS(2,:),re_GPS(3,:));
% end
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');
% scatter3(0,0,0,250,'b','filled');
% legend(GPS.textdata(2:end),'Earth ');
% title('Orbit of GPS in ECEF coordinate system')

figure
view (3)
hold on
grid on
Orbit=[];
for i=1:length(GLONASS.data)
    Orbit=Position_GLONASS(:,(i-1)*size(t,2)+1:i*size(t,2));
    plot3(Orbit(1,:),Orbit(2,:),Orbit(3,:));
end
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
scatter3(0,0,0,250,'b','filled');

title('Orbit of GLONASS in geocentric coordinate system')


% figure
% view (3)
% hold on
% grid on
% for i=1:length(GLONASS.data)
%     Orbit=Position_GLONASS(:,(i-1)*size(t,2)+1:i*size(t,2));
%     re_GLONASS =cell2mat( arrayfun(@(p) rot_GLONASS(:,:,p)*Orbit(:,p), 1:length(t), 'UniformOutput', false));
%     plot3(re_GLONASS(1,:),re_GLONASS(2,:),re_GLONASS(3,:));
% end
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');
% scatter3(0,0,0,250,'b','filled');
% title('Orbit of GLONASS in ECEF coordinate system')

figure
view (3)
hold on
grid on
Orbit=[];
for i=1:length(Galileo.data)
    Orbit=Position_Galileo(:,(i-1)*size(t,2)+1:i*size(t,2));
    plot3(Orbit(1,:),Orbit(2,:),Orbit(3,:));
end
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
scatter3(0,0,0,250,'b','filled');

title('Orbit of Galileo in geocentric coordinate system')

% figure
% hold on% view (3)

% grid on
% for i=1:length(Galileo.data)
%     Orbit=Position_Galileo(:,(i-1)*size(t,2)+1:i*size(t,2));
%     re_Galileo =cell2mat( arrayfun(@(p) rot_Galileo(:,:,p)*Orbit(:,p), 1:length(t), 'UniformOutput', false));
%     plot3(re_Galileo(1,:),re_Galileo(2,:),re_Galileo(3,:));
% end
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');
% scatter3(0,0,0,250,'b','filled');
% title('Orbit of Galileo in ECEF coordinate system')