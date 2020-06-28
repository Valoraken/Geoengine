% Physical Gerodesy Lab 3 
% Xiao, Tianqi 3371477
% Wang, Zhenqiao 3371590
% 30.05.2019
% last edit: 6.6.2019

close all
clear
clc

% Constant and corrdinates
k=77 + 90                                 % k value 
r = 6371*1e3;                             % spherical coordinate r [m] 
phi=[deg2rad(20+2*k) deg2rad(0:0.5:90)];  % spherical coordinate phi [rad]
rho= 5515;                                % density of the earth [kg/m^3]
lamda=deg2rad(10);                        % spherical coordinate lamda [rad]
w= 7.292115*1e-5;                         % angular velocity of the Earth [s^-1]
G= 6.672*1e-11;                           % gravitational constant
RE=6371*1e3;                              % radius of the Earth [m]
x=[r.*cos(phi).*cos(lamda);r.*cos(phi).*sin(lamda);r.*sin(phi)];    % Cartesian coordinates of the point
%% Task 1
% Gravitational Potential
V=4/3*pi*G*rho*RE.^3/r;
Vc=1/2*w^2*r.^2.*(cos(phi)).^2;
W=V+Vc;
% Attraction
a=-4/3/r^3*pi*G*rho*RE^3*x;
ac=w^2*[x(1:2,:);zeros(1,length(x))];
g=a+ac;

xi=(acos(dot(a,g)./(vecnorm(a).*vecnorm(g))));% [rad]
delta_g=vecnorm(g)-vecnorm(a);

fprintf('Numerical results:\n ');
fprintf('Gravitational potential:\nV = %.6f m^2/s^2 \n ',V(1));
fprintf('Centrifugal potential:\nVc = %.6f m^2/s^2 \n ',Vc(1));
fprintf('Gravity potential:\nW = %.6f m^2/s^2 \n ',W(1));
fprintf('Gravitational attraction:\n   [%.4f\na = %.4f  m/s^2\n     %.4f]\n',a(:,1));
fprintf('Centrifugal acceleration:\t\n    [%.4f\nac = %.4f  m/s^2\n     %.4f]\n',ac(:,1));
fprintf('Gravity attraction:\ng = %.4f m/s^2 \n ',norm(g(:,1)));
fprintf('Disturbance of the direction:\n¦Î = %.6f¡ã\n ',rad2deg(xi(1)));
fprintf('Disturbance of the attraction:\ndg = %.6f m/s^2 \n ',delta_g(1));

% plot
figure(1)
hold on
plot(rad2deg(phi(2:length(phi))),delta_g(2:length(phi)))
title('Variation Tendency of Attraction (Absolute Value)');
xlabel('Latitude \Phi [degree]');
ylabel('Absolute value of attraction \delta g [m/s^{2}]');

figure(2)
hold on
plot(rad2deg(phi(2:length(phi))),Vc(2:length(phi)))
title('Variation Tendency of Centrifugal Potential');
xlabel('Latitude \Phi [degree]');
ylabel('Centrifugal potential Vc [m^{2}/s^{2}]');

figure(3)
hold on
plot(rad2deg(phi(2:length(phi))),rad2deg(xi(2:length(phi))))
title('Variation Tendency of Direction Disturbances');
xlabel('Latitude \Phi [degree]');
ylabel('Direction \xi [degree]');
%% Task 2
fprintf('\nTask 2\n ')
Phi= 42;                   % latitude [deg]
theta= deg2rad(90-Phi);     % co-latitude [rad]
v= 400*1000/3600;            % velocity [m/s]
% Coriolis acceleration
% in this task vE = vN = v
vE= v; vN= v;
a_EW= 2*w*[-cos(theta)*vE;0;sin(theta)*vE];
a_NS= 2*w*[0;cos(theta)*vN;0];

fprintf('Coriolis acceleration for EW direction:\n     [%.4f\na_EW = %.4f  m/s^2\n       %.4f]\n',a_EW);
fprintf('Coriolis acceleration for NS direction:\n      [%.4f\na_NS = %.4f  m/s^2\n       %.4f]\n\n',a_NS);
% Eoetvoes correction
E_EW=-2*w*sin(theta)*vE;
E_NS=0;
fprintf('Eoetvoes correction for EW direction:\nE_EW = %.4f m/s^2\n',E_EW);
fprintf('Eoetvoes correction for NS direction:\nE_NS = %.4f m/s^2\n\n',E_NS);

% Calculate accuracy via error propagation
sigma_E=1e-5;
sigma_v=sqrt(sigma_E^2/(2*w*sin(theta))^2);
fprintf('Accuracy of velocity:\nSigma_v = %.4f m/s\n',sigma_v);
