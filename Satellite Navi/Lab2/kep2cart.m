function [Position,Velocity] = kep2cart(a,e,I,omega,Omega,M_0,T0,t)
% Function kep2cart coverts Keplerian elements (K.E.) to Cartesian coordinates and velocities
% Computation of the state vector (position and velocity) of a satellite in a given epoch, from the six keplerian elements. The state vector is refered to the Earth Centered Inertial (ECI) Reference Frame.
% Input:        
%         a             [1 x 1] semi-major axis in meters.
%         e             [1 x 1] eccentricity.
%         I             [1 x 1] inclination in rad.
%         omega         [1 x 1] argument of perigee in rad.
%         Omega         [1 x 1] longitude of ascending node in rad.
%         T0            [1 x 1] epoch of reference in second.
%         t             [n x 1] epoch of computation in second.
%
% Output: Positiom      [3 x n] [X Y Z]'coordinates  in meters.
%         Velocity      [3 x n] [VX VY VZ]' velocity in m/s.
%


%% Define GM
GM = 3.98600415*1e14;
%% Calculate parameters
% mean motion
n = sqrt(GM./a.^3);
% mean anomaly
M = M_0+n.*(t-T0);
%Normalize mean anomaly
while max(M >= 2*pi)
    M(M >= 2*pi) = M(M >= 2*pi)-2*pi; 
end
while max(M < 0)    
    M(M < 0)     = M(M < 0)    +2*pi; 
end

% iteration to get eccentric anomaly
E0=M; % E0
E= M+e.*sin(E0);
Ei=E0;
while max(abs(Ei-E))>1e-14
    Ei = E;
    E = M+e.*sin(Ei);
end


%% Get Position and Velocity in the orbital plane coordinate system
% Position 
x = a.*(cos(E)-e);
y = a.*sqrt(1-e.^2).*sin(E);
z = zeros(size(x));

% Velocity
r = a.*(1-e.*cos(E));
vx = -n.*(a.^2./r).*sin(E);
vy = n.*(a.^2./r).*sqrt(1-e.^2).*cos(E);
vz = zeros(size(vx));

%% Transfer into geocentric coordinate system:
% Position = angle2dcm( -omega, -I, -Omega, 'ZXZ' )*[x;y;z];
% Velocity = angle2dcm( -omega, -I, -Omega, 'ZXZ' )*[vx;vy;vz];
Position = cell2mat( arrayfun(@(p) angle2dcm( -omega(p), -I(p), -Omega(p), 'ZXZ' )*[x(p,:);y(p,:);z(p,:)], 1:size(x,1), 'UniformOutput', false));
Velocity = cell2mat( arrayfun(@(p) angle2dcm( -omega(p), -I(p), -Omega(p), 'ZXZ' )*[vx(p,:);vy(p,:);vz(p,:)], 1:size(x,1), 'UniformOutput', false));

end

