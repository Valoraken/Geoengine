% Physical Geodesy Exercise 2
% NAME: Xiao, Tianqi, Wang Zhenqiao
% Student ID: 3371477, 3371590
% 25.5.2019

clear 
clc
close all

% K-value
k = (90+77)/2

%% Set Numerical Values
mE = 5.9736e24;               % mass of the Earth [kg]
RE = 6371000;                 % radius of the Earth [m]
XE = 0;YE= 0;                 % coordinates of the Earth wrt. a geocentric reference system [m]
mM = 7.349e22;                % mass of the Moon [kg]
RM = 1738000;                 % radius of the Earth [m]
rM = 384400000;               % Distance Earth-Moon [m]
XM = rM*cos(deg2rad(k*10));   % coordinates of the moon wrt. a geocentric reference system [m]
YM = rM*sin(deg2rad(k*10));
G=6.672e-11;                  % gravitational constant[m3s-2kg-1]

% Calculate density
rhoE = mE/(4/3*pi*RE^3);      % density of the Earth
rhoM = mM/(4/3*pi*RM^3);      % density of the Moon
%% Part 1: Superposition of Earth's and Moon's gravitational fields

% build dense grid
[X Y] = meshgrid(linspace(-rM,0.5*rM,800),linspace(-0.5*rM,rM,800));

% superposition of potentials
VE = V_sphere(X,Y,RE,XE,YE,rhoE);
VM = V_sphere(X,Y,RM,XM,YM,rhoM);
V = VE+VM;
V = reshape(V,size(X));

%visualize Potential
figure (1)
hold on
[CV,hV]=contour(X,Y,V/1e5,[2:1:16,16:2:24,24:7:80,80:50:180,180:80:500],'b ');
clabel(CV,hV,'FontSize',8,'Color','red');
axis equal;
xlabel('x [m]');
ylabel('y [m]');
title('Potential of both Earth and Moon (*e^5)');

% attraction
[X Y] = meshgrid(linspace(-rM,0.5*rM,15),linspace(-0.5*rM,rM,15));
[axE ayE] = a_sphere(X,Y,RE,XE,YE,rhoE);
[axM ayM] = a_sphere(X,Y,RM,XM,YM,rhoM);
ax = axE+axM;
ay = ayE+ayM;
ax = reshape(ax,size(X));
ay = reshape(ay,size(X));
% normalize
amax = max(max(abs(ax(:))),max(abs(ay(:))));
ax = ax/amax;
ay = ay/amax;

% visualize
figure (2)
hold on
scatter(XM,YM,'r');
text(XM+1e7,YM+1e7,'Moon');
scatter(XE,YE,'b');
text(XE+1e7,YE+3e7,'Earth');
quiver(X,Y,ax,ay,1.2);
grid on
axis equal
xlabel('x [m]');
ylabel('y [m]');
title('Attraction of both Earth and Moon (Vector Field)');

%% Part 2 : Gravitational potential and attraction of spherical shells

% Numerical Values
Rc = 3500000;                % core radius[m]
rho_c = 11200;               % core density [kg/m3]
Rm = 6400000;                % mantle radius[m]
rho_m = 4300;                % mantle density[kg/m3]

r = 0:100000:4*Rm;
% for shell
Vshell = V_shell(Rc,Rm,rho_m,r);
ashell = a_shell(Rc,Rm,rho_m,r);
% for sphere
Vsphere = V_sphere(0,r,Rc,0,0,rho_c);
[a_sx,a_sy] = a_sphere(0,r,Rc,0,0,rho_c);
asphere = -sqrt(a_sx.^2+a_sy.^2);

figure(3)
subplot(2,1,1)
hold on
plot(r,Vshell+Vsphere,'b');
xlabel('distance [m]');
ylabel('potential [J]');
title('Potential of Simplified Model of Earth');

subplot(2,1,2)
plot(r,ashell+asphere,'m ');
xlabel('distance [m]');
ylabel('attraction [m/s^2]');
title('Attraction of Simplified Model of Earth');

%% Part 3 PREM density model of the Earth

% load and plot density distribution
load('PREM.mat')
figure(4)
plot(PREM(:,1)*1000,PREM(:,2));
xlabel('Radial coordinate r [m]');
ylabel('Density ¦Ñ [kg/m^3]');
title('PREM');

r = 0:10:2*RE/1000;
% The core can be seen as a sphere
V0 = V_sphere(0,r,10,0,0,PREM(1,2));
[a0x,a0y] = a_sphere(0,r,10,0,0,PREM(1,2));
a0 = -sqrt(a0x.^2+a0y.^2);

% Then every 10km-layer can be seen as a shell
i=2;
for r_in = 10:10:RE/1000-10
     r_out = r_in+10;
     % find a midium value of density in PREM
     rho = PREM(find(ceil(PREM(:,1))==ceil((r_in+r_out)/2)),2);
     % Calculate potential and attraction
     Vs(i,:) = V_shell(r_in,r_out,rho,r);
     as(i,:) = a_shell(r_in,r_out,rho,r);
     i = i+1;
end
 Vs(1,:) = V0;
 V_PREM = sum(Vs);
 as(1,:) = a0;
 a_PREM = sum(as);
 
% visualize
figure
subplot(2,1,1)
hold on
plot(r,V_PREM*10^6)
scatter(find(max(V_PREM))*10-10,max(V_PREM)*(10^6),'p');
legend('potential','the largest potential');
xlabel('distance [km]');
ylabel('potential [J]');
title('potential of earth (PREM)');

subplot(2,1,2)
hold on
[maxa,maxa_index] = max(abs(a_PREM));
plot(r,a_PREM*10^3,'m')
scatter(maxa_index*10-10,a_PREM(maxa_index)*1000,'p');
legend('attraction','the largest attraction');
xlabel('distance [km]');
ylabel('attraction [m/s^2]');
title('attraction of earth (PREM)');

% Print results
fprintf('Task 9: \n\n');
fprintf('the gravitaional attraction at r = RE is %.4f [m/s^2] \n',a_PREM(ceil(RE/10000))*1000);
fprintf('the gravitaional potential at r = RE is %.6s [J]\n\n',V_PREM(ceil(RE/10000))*10^6);

fprintf('Task 10: \n\n');
fprintf('the largest gravitaional attraction is %.4f [m/s^2] \n',maxa*1000);
fprintf('the position of the largest gravitaional attraction is %d km\n',maxa_index*10);
fprintf('the corresponding gravitaional potential is %.6s [J]\n\n',V_PREM(maxa_index)*10^6);
fprintf('the largest gravitaional potential is %.6s [J] \n',max(V_PREM)*10^6);
fprintf('the position of the largest gravitaional potential is %d [m]\n',find(max(V_PREM))*10-10);
fprintf('the corresponding gravitaional attraction is %.4f [m/s^2]\n\n',a_PREM(find(max(V_PREM)))*1000);