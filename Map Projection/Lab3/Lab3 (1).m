% Xiao, Tianqi 3371477
% Lab3 16.12.18
% last changes 18.12.18

% This program is to realize transformation in the Meta-Coorinate System
% and compute meta-coordinates of the map for the locations: Erfurt, Tokyo
% and Washington D.C  with respected to both meta poles
% Meta pole: Beijing(116��23'E ,39��56'N); Perth(115��51'E ,31��57'N)
close all
clear all
clc
%% load data
load coast.mat
% from degree to [rad]
lamda=long*pi/180;
phi=lat*pi/180;
% set R [m] and phi0 [deg] for Stereographic projection
R=1;
phi0=60;

% spherical corrinates:
% Beijing [rad]
lamda01=deg2rad(dms2degrees([116 23 0]));
phi01=deg2rad(dms2degrees([39 56 0]));
% Perth [rad]
lamda02=deg2rad(dms2degrees([115 51 0]));
phi02=-deg2rad(dms2degrees([31 57 0]));
% Erfurt Tokyo and Washington D.C. [rad]
Lamda_E=deg2rad(dms2degrees([11 02 0]));     %E
Lamda_T=deg2rad(dms2degrees([139 46 0]));    %E
Lamda_W=-deg2rad(dms2degrees([77 01 0]));    %W

Phi_E=deg2rad(dms2degrees([50 58 0]));   %N
Phi_T=deg2rad(dms2degrees([35 42 0]));   %N
Phi_W=deg2rad(dms2degrees([38 54 0]));   %N
%% example data in seminar 5
% lamda0=deg2rad(dms2degrees([9 11 0]));
% phi0=deg2rad(dms2degrees([48 46 0]))
% omega0=deg2rad(30.75)

%% find omega0
% After the rotation, Stuttgart is lying on the meta-Greenwich meridian
% the meta coordinates of Stuttgart (A_s,B_s) 
% with A_s=arctan(sin(��lamda)/(sin(phi0)*cos(��lamda)-cos(phi0)*tan(phi))-omega0=0
% omega0=arctan(sin(��lamda)/(sin(phi0)*cos(��lamda)-cos(phi0)*tan(phi)

% Stuttgart [rad]
lamda_s=deg2rad(dms2degrees([9 11 0]));
phi_s=deg2rad(dms2degrees([48 46 0]));

% omega0 [rad]
omega01=atan2(sin(lamda_s-lamda01),(sin(phi01)*cos(lamda_s-lamda01)-cos(phi01)*tan(phi_s)))
omega02=atan2(sin(lamda_s-lamda02),(sin(phi02)*cos(lamda_s-lamda02)-cos(phi02)*tan(phi_s)))

%% Rotation method for the pole Beijing
% rotation matrix
Romega01=[cos(omega01) sin(omega01) 0;-sin(omega01) cos(omega01) 0;0 0 1];
Rdelta01=[cos(pi/2-phi01) 0 -sin(pi/2-phi01);0 1 0 ;sin(pi/2-phi01) 0 cos(pi/2-phi01)];
Rlamda01=[cos(lamda01) sin(lamda01) 0;-sin(lamda01) cos(lamda01) 0;0 0 1];

% Cartesian coordinates before rotation [m]
X_C=[cos(lamda).*cos(phi) sin(lamda).*cos(phi) sin(phi)]';
% Cartesian coordinates after rotation [m]
X_M=(Romega01*Rdelta01*Rlamda01*X_C)';

% points in meta-coorinates [rad]
A_Beijing=atan2(X_M(:,2),X_M(:,1));
B_Beijing=asin(X_M(:,3));

%% direct formula for the pole Perth
% points in meta-coorinates [rad]
B_Perth=asin(sin(phi).*sin(phi02)+cos(phi).*cos(phi02).*cos(lamda-lamda02));
A_Perth=atan2(sin(lamda-lamda02),(sin(phi02).*cos(lamda-lamda02)-cos(phi02).*tan(phi)))-omega02;

%% Computation of meta-coorinates 
Lamda=[Lamda_E;Lamda_T;Lamda_W];
Phi=[Phi_E;Phi_T;Phi_W];
% Rotaion method for pole Beijing
C=[cos(Lamda).*cos(Phi) sin(Lamda).*cos(Phi) sin(Phi)]';
M=(Romega01*Rdelta01*Rlamda01*C)';
DegA_Beijing=rad2deg(atan2(M(:,2),M(:,1)));
DegB_Beijing=rad2deg(asin(M(:,3)));
DmsA_Beijing=degrees2dms(DegA_Beijing);
DmsB_Beijing=degrees2dms(DegB_Beijing);

fprintf('With respect to pole Beijing the meta-coordinates (A ,B) for the locations are: \nin degrees:\n');
fprintf('Erfurt: %.4f E   %.4f N\n',DegA_Beijing(1),DegB_Beijing(1));
fprintf('Tokyo: %.4f E   %.4f N\n',DegA_Beijing(2),DegB_Beijing(2));
fprintf('Washington D.C.: %.4f E   %.4f N\n',DegA_Beijing(3),DegB_Beijing(3));
fprintf('\nin dms: \n');
fprintf('Erfurt: %d %d %.4f E   %d %d %.4f N\n',DmsA_Beijing(1,:),DmsB_Beijing(1,:));
fprintf('Tokyo: %d %d %.4f E   %d %d %.4f N\n',DmsA_Beijing(2,:),DmsB_Beijing(2,:));
fprintf('Washington D.C.: %d %d %.4f E   %d %d %.4f N\n',DmsA_Beijing(3,:),DmsB_Beijing(3,:));

% Direct formula for pol Perth
DegB_Perth=rad2deg(asin(sin(Phi).*sin(phi02)+cos(Phi).*cos(phi02).*cos(Lamda-lamda02)));
DegA_Perth=rad2deg(atan2(sin(Lamda-lamda02),(sin(phi02).*cos(Lamda-lamda02)-cos(phi02).*tan(Phi)))-omega02);
DegA_Perth(DegA_Perth>180)=DegA_Perth(DegA_Perth>180)-360;
DmsA_Perth=degrees2dms(DegA_Perth);
DmsB_Perth=degrees2dms(DegB_Perth);

fprintf('\nWith respect to pole Perth the meta-coordinates (A ,B) for the locations are: \nin degrees:\n');
fprintf('Erfurt: %.4f E   %.4f N\n',DegA_Perth(1),DegB_Perth(1));
fprintf('Tokyo: %.4f E   %.4f N\n',DegA_Perth(2),DegB_Perth(2));
fprintf('Washington D.C.: %.4f E   %.4f N\n',DegA_Perth(3),DegB_Perth(3));
fprintf('\nin dms: \n');
fprintf('Erfurt: %d %d %.4f E   %d %d %.4f N\n',DmsA_Perth(1,:),DmsB_Perth(1,:));
fprintf('Tokyo: %d %d %.4f E   %d %d %.4f N\n',DmsA_Perth(2,:),DmsB_Perth(2,:));
fprintf('Washington D.C.: %d %d %.4f E   %d %d %.4f N\n',DmsA_Perth(3,:),DmsB_Perth(3,:));
%% Computation of the map (Cartesian) 
% pole Beijing [m]
[B_x,B_y]=Stereographic(DegA_Beijing,DegB_Beijing,phi0,R); 
fprintf('\nWith respect to pole Beijing the Cartesian coordinates of the map(x,y) for the locations are: \nin degrees:\n');
fprintf('Erfurt: %.3f ,%.3f [m]\n',B_x(1),B_y(1));
fprintf('Tokyo: %.3f ,%.3f [m]\n',B_x(2),B_y(2));
fprintf('Washington D.C.: %.3f ,%.3f [m] \n',B_x(3),B_y(3));

% pole Perth [m]
[P_x,P_y]=Stereographic(DegA_Perth,DegB_Perth,phi0,R); 
fprintf('\nWith respect to pole Perth the Cartesian coordinates of the map(x,y) for the locations are: \nin degrees:\n');
fprintf('Erfurt: %.3f ,%.3f [m]\n',P_x(1),P_y(1));
fprintf('Tokyo: %.3f ,%.3f [m]\n',P_x(2),P_y(2));
fprintf('Washington D.C.: %.3f ,%.3f [m] \n',P_x(3),P_y(3));

%% virtualize
% from rad to degree
long_B=rad2deg(A_Beijing);
lat_B=rad2deg(B_Beijing);

long_P=rad2deg(A_Perth);
lat_P=rad2deg(B_Perth);

% avoid phi<0, and lamda>360
[long_B, lat_B] = jump2nan(long_B, 180,lat_B);
lat_B(lat_B<0)=NaN;

long_P(long_P>180) = long_P(long_P>180)-360;
[long_P, lat_P] = jump2nan(long_P, 180,lat_P);
lat_P(lat_P<0)=NaN;

% Stereographic projection
[x_B,y_B]=Stereographic(long_B,lat_B,phi0,R); 
[x_P,y_P]=Stereographic(long_P,lat_P,phi0,R); 

% the parameter lines within the interval
% the interval of parameter lines
long_interval=-180:30:180;
lat_interval=-90:15:90;
% meridians 
[long_m,lat_m]=meshgrid(long_interval ,linspace(min(lat_interval),max(lat_interval),360)); 
long_m=deg2rad(long_m);
lat_m=deg2rad(lat_m);
B_mB=rad2deg(asin(sin(lat_m).*sin(phi01)+cos(lat_m).*cos(phi01).*cos(long_m-lamda01)));
A_mB=rad2deg(atan2(sin(long_m-lamda01),(sin(phi01).*cos(long_m-lamda01)-cos(phi01).*tan(lat_m)))-omega01);
B_mP=rad2deg(asin(sin(lat_m).*sin(phi02)+cos(lat_m).*cos(phi02).*cos(long_m-lamda02)));
A_mP=rad2deg(atan2(sin(long_m-lamda02),(sin(phi02).*cos(long_m-lamda02)-cos(phi02).*tan(lat_m)))-omega02);
% [A_mB, B_mB] = jump2nan(A_mB, 180,B_mB);
B_mB(B_mB<0)=NaN;

A_mP(A_mP>180) = A_mP(A_mP>180)-360;
% [A_mP, B_mP] = jump2nan(A_mP, 180,B_mP);
B_mP(B_mP<0)=NaN;
[x_mB,y_mB]=Stereographic(A_mB,B_mB,phi0,R);
[x_mP,y_mP]=Stereographic(A_mP,B_mP,phi0,R);

% parallels
[long_p,lat_p]=meshgrid(linspace(min(long_interval),max(long_interval),360),lat_interval);

long_p=deg2rad(long_p);
lat_p=deg2rad(lat_p);
B_pB=rad2deg(asin(sin(lat_p).*sin(phi01)+cos(lat_p).*cos(phi01).*cos(long_p-lamda01)));
A_pB=rad2deg(atan2(sin(long_p-lamda01),(sin(phi01).*cos(long_p-lamda01)-cos(phi01).*tan(lat_p)))-omega01);
B_pP=rad2deg(asin(sin(lat_p).*sin(phi02)+cos(lat_p).*cos(phi02).*cos(long_p-lamda02)));
A_pP=rad2deg(atan2(sin(long_p-lamda02),(sin(phi02).*cos(long_p-lamda02)-cos(phi02).*tan(lat_p)))-omega02);
% [A_mB, B_mB] = jump2nan(A_mB, 180,B_mB);
B_pB(B_pB<0)=NaN;

A_pP(A_pP>180) = A_pP(A_pP>180)-360;
% [A_mP, B_mP] = jump2nan(A_mP, 180,B_mP);
B_pP(B_pP<0)=NaN;
[x_pB,y_pB]=Stereographic(A_pB,B_pB,phi0,R);
[x_pP,y_pP]=Stereographic(A_pP,B_pP,phi0,R);

%edges
[ex,ey]=Stereographic(-180:180,0,phi0,R);
% desplay the results and save
% meta pole Beijing
figure(1)
hold on
plot(x_B,y_B)
plot(x_mB,y_mB,'color', [0.5 0.5 0.5]);
plot(x_pB',y_pB','color', [0.5 0.5 0.5]);
plot(ex',ey','color', [0.5 0.5 0.5]);
scatter(0,0,8,'b');
scatter(B_x(1:2),B_y(1:2),10,'r');
text(0,0,'Beijing','Color','b','fontsize',14)
text(B_x(1),B_y(1),' Erfurt','Color','r','fontsize',14)
text(B_x(2),B_y(2),' Tokyo','Color','r','fontsize',14)
title('Sterographic Projection with Meta Pole Beijing ');
axis equal
axis off
print('-depsc2', 'Beijing.eps');
print('-dpng', '-r800', 'Beijing.png')

% meta pole Beijing
figure(2)
hold on
plot(x_P,y_P)
plot(x_mP,y_mP,'color', [0.5 0.5 0.5]);
plot(x_pP',y_pP','color', [0.5 0.5 0.5]);

scatter(0,0,8,'b');
scatter(P_x(2),P_y(2),10,'r');
text(0,0,' Perth','Color','b','fontsize',14)
text(P_x(2),P_y(2),' Tokyo','Color','r','fontsize',14)
title('Sterographic Projection with Meta Pole Perth ');
axis equal
axis off
print('-depsc2', 'Perth.eps');
print('-dpng', '-r800', 'Perth.png')

