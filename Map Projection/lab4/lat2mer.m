function G = lat2mer(A, E2, B)
% Meridional arc length G from equator until ellipsoidal latitude B
% IN:
%    A ..... semi major axis of reference ellipsoid in [m]
%    E2 .... square of 1st numerical excentricity of 
%            reference ellipsoid in [-]
%    B ..... ellipsoidal latitude in [rad]
% OUT:
%    G ..... meridional arc length in [m]

% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GI, Uni Stuttgart
% -------------------------------------------------------------------------
% revision history:
%    2016-01-19: MR, brush up help text
%    2011-12-20: MR, inital version after "geomed_merbod.pdf"
% -------------------------------------------------------------------------


e_0  = 1 - 1/4 * E2 - 3/64 * E2^2 -   5/256 * E2^3 - 175/16384 * E2^4 -   441/65536 * E2^5;
e_2  =    -3/8 * E2 - 3/32 * E2^2 - 45/1024 * E2^3 -  105/4096 * E2^4 - 2205/131072 * E2^5;
e_4  =              15/256 * E2^2 + 45/1024 * E2^3 + 525/16384 * E2^4 +  1575/65536 * E2^5;
e_6  =                             -35/3072 * E2^3 - 175/12288 * E2^4 - 3675/262144 * E2^5;
e_8  =                                              315/131072 * E2^4 + 2205/524288 * E2^5;
e_10 =                                                                 -693/1310720 * E2^5;

G   = A * (e_0 * B + e_2 * sin(2 * B) + e_4 * sin(4 * B) + e_6 * sin(6 * B) + e_8 * sin(8 * B) + e_10 * sin(10 * B));
