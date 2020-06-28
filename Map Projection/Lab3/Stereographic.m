function [x,y] = Stereographic(longtitude,latitude,phi0,R)
 % This function calculates the point[x,y] in the chart of Stereographic projection
 % Gnomonic Projection
 % Input:
 % longtitude,latitude are the spherical coordinates (the unit is degree)
 % phi0 = const, 0 <= phi0 <= 90 (the unit is degree)
 %
 % R is the radius of sphere. (the unit is kilometer)
 % Output:
 % [x,y] is the coordinates in the plane.
 
 % Check input

if nargin < 4
  error('Not enough arguments');
end
if  (phi0>90 || phi0<0)
    error('phi0 should be a const between 0~90');
end

 % Change the degree into rad, caculate Lamda and Delta
 Lamda = longtitude * pi / 180;
 Delta = (90-latitude) * pi /180;
 
 % The coordinates in the plane
 x = R*(1+sin(phi0*pi/180))*tan(Delta/2).*cos(Lamda);
 y = R*(1+sin(phi0*pi/180))*tan(Delta/2).*sin(Lamda);
    
end


