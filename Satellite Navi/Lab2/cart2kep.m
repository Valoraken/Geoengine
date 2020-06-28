function [a,e,I,W,w] = cart2kep(Position , Velocity)
% Function cart2kep coverts Cartesian coordinates and velocities to Keplerian elements (K.E.)
% Computation of the five keplerian elements form the state vector (position and velocity) of a satellite refered to the Earth Centered Inertial (ECI) Reference Frame in a given epoch
% Input:        
%         Positiom      [3 x n] [X Y Z]'coordinates  in meters.
%         Velocity      [3 x n] [VX VY VZ]' velocity in m/s.

%
% Output: 
%         a              semi-major axis in meters.
%         e              eccentricity.
%         I              inclination in rad.
%         w              argument of perigee in rad.
%         W              longitude of ascending node in rad.

%% Define GM
GM = 3.98600415*1e14;
%%
h=cross(Position,Velocity);
W=atan2(h(1,:),-h(2,:));
I=atan2(sqrt(h(1,:).^2+h(2,:).^2),h(3,:));
r = sqrt(Position(1,:).^2+Position(2,:).^2+Position(3,:).^2);
v = sqrt(Velocity(1,:).^2+Velocity(2,:).^2+Velocity(3,:).^2);
a = 1./((2./r)-v.^2./GM);
n = sqrt(GM./a.^3);
e = sqrt(1-(h(1,:).^2+h(2,:).^2+h(3,:).^2)./(a.*GM));
E = atan2(dot(Position,Velocity)./(a.^2.*n),(1-r./a));
niu = 2*atan2(sqrt((1+e)./(1-e)),1./tan(E./2));
u = acos((Position(1,:).*cos(W)+Position(2,:).*sin(W))./r);
u(find(Position(3,:)<0))=2*pi-u(find(Position(3,:)<0));
w = u-niu;
w(w<0)=w(w<0)+2*pi;
% plot(rad2deg(w));

end

