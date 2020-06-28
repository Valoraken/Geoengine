% Physical Geodesy Exercise 2
% NAME: Xiao, Tianqi, Wang Zhenqiao
% Student ID: 3371477, 3371590
% 23.5.2019

function [ax,ay] = a_sphere(X,Y,R,Xm,Ym,rho)

% This function is to calculate the vector valued gravitational attractions of a
% homogeneous solid sphere at a given certain point
%    INPUT  :
% ( X,Y ) : coordinates of the given point
%    R    : radius of the solid sphere
% (Xm,Ym) : coorinates of the mass centre of the sphere
%   rho   : density of the solid sphere
%    OUTPUT :
% (ax,ay) : gravitational attracion at the point (X,Y)

%% Define constant and calculate distance
G= 6.672e-11;                   % gravitational constant
r= sqrt((X-Xm).^2+(Y-Ym).^2);     % distance between mass centre and given point

%% Calculate Attraction
a(r<=R)= -4/3*pi*G.*rho.*r(r<=R);
a(r>R)= -4/3*pi*G.*rho*R^3./r(r>R).^2;
alpha = atan2(Y-Ym,X-Xm);
ax =cos(alpha(:)').*a;
ay = sin(alpha(:)').*a;
end

