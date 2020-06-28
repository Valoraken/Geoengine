% Physical Geodesy Lab 3
% Task 2 and Task 3
% Name : Yi Hong/Shu Suo/Chanjuan Zhang
% Matriculation number: 3294211/3295472/3296099
%
clear
clc
%% Task 2 
k = 9;                                                     % k value 
G = 6.672*10^(-11);                                        % gravitational constant [m^3*s^(-2)*kg^(-1)]
R = 6371*10^3;                                             % radium of earth [m] 
Rho = 5515;                                                % [kg*m^(-3)]
r = R;
omega = 7.292115*10^(-5);

% Quantities at point p
lambda_p = 10;
phi_p = 20+2*k;
thet_p = 90-phi_p;   
x_p = r*sin(deg2rad(thet_p))*cos(deg2rad(lambda_p));
y_p = r*sin(deg2rad(thet_p))*sin(deg2rad(lambda_p));
z_p = r*cos(deg2rad(thet_p));

V_p = 4*pi*G*Rho*R^3/(3*r);
Vc_p = omega^2*r^2*sin(deg2rad(thet_p))^2/2;
W_p = Vc_p + V_p;

a = -4*pi*G*Rho*r/3;                                        % magnitude of a
a_x = a*x_p/r;
a_y = a*y_p/r;
a_z = a*z_p/r;
a_p = [a_x; a_y; a_z];
ac_p = omega^2.*[x_p;y_p;0];
g_p = ac_p+a_p;
g_p_norm = norm(g_p);                                       % magnitude of g
 
delta_W = W_p - V_p;                                        % disturbance of the potential
delta_g_p = g_p_norm - norm(a_p);                           % absolute value of the attraction
xi_p = asin(norm(ac_p)*sin(deg2rad(phi_p))/g_p_norm);       % disturbance of the direction (sine law) [rad]
xi_p = rad2deg(xi_p);                                       % [deg]

% Sketches
phi = [0:5:90];
thet = 90-phi;
lambda = 20;
x = r*sin(deg2rad(thet))*cos(deg2rad(lambda));
y = r*sin(deg2rad(thet))*sin(deg2rad(lambda));
z = r*cos(deg2rad(thet));

figure(1);
Vc = omega^2*r^2*sin(deg2rad(thet)).^2/2;
plot(phi,Vc);
title('Variation Tendency of Centrifugal Potential');
xlabel('Latitude \Phi [degree]');
ylabel('Centrifugal potential Vc [m^{2}/s^{2}]  ');

figure(2);
a = -4*pi*G*Rho*r/3;                              
a_x = a*x/r;
a_y = a*y/r;
a_z = a*z/r;
a = [a_x; a_y; a_z];
ac = omega^2.*[x;y;zeros(size(x))];
g = ac + a;
[m, n] = size(g);
for i = 1:n
    delta_g(i) = norm(g(:,(i))) - norm(a(:,(i)));
end
plot(phi,delta_g);
title('Variation Tendency of Attraction (Absolute Value)');
xlabel('Latitude \Phi [degree]');
ylabel('Absolute value of attraction \delta g [m/s^{2}]');

figure(3);
for i = 1:n
    xi(i) = asin(norm(ac(:,(i)))*sin(deg2rad(phi(i)))/norm(g(:,(i)))); 
end
xi = rad2deg(xi);
plot(phi,xi);
title('Variation Tendency of Direction Disturbances');
xlabel('Latitude \Phi [degree]');
ylabel('Direction \xi [degree]');

%% Task 3
% constant
phi=42/180*pi;
theta=pi/2-phi;
v=400*1000/3600;
omega=7.292115e-5;
v_EW=[-v,0]';                                                 % velocity of east to west, posotive direction is east
v_NS=[0,-v]';                                                 % velocity of north to south, positive direction is north

% coriolis acceleration
a_cor_EW_t=2*omega*[-cos(theta), cos(theta)]*v_EW;
a_cor_NS_t=2*omega*[-cos(theta), cos(theta)]*v_NS;
% Eötvös correction
eoet_correction_EW=-2*omega*sin(theta)*v_EW(1,1);
eoet_correction_NS=-2*omega*sin(theta)*v_NS(1,1);

% accuraycy of velocity
sigma_eoet=1e-5; % 1mGal=10^-5 m/s^2
sigma_v=sigma_eoet/(2*omega*sin(theta));