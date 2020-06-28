function B = mer2lat(A, E2, G)
% Ellipsoidal latitude B for given meridional arc length G
% IN:
%    A ..... semi major axis of reference ellipsoid in [m]
%    E2 .... square of 1st numerical excentricity of 
%            reference ellipsoid in [-]
%    G ..... meridional arc length in [m]
% OUT:
%    B ..... ellipsoidal latitude in [rad]

% -------------------------------------------------------------------------
% authors:
%    Matthias ROTH (MR), GI, Uni Stuttgart
% -------------------------------------------------------------------------
% revision history:
%    2016-01-19: MR, brush up help text
%    2011-12-20: MR, inital version after "geomed_merbod.pdf"
% -------------------------------------------------------------------------

e_0  = A * (1 - 1/4 * E2 - 3/64 * E2^2 -    5/256 * E2^3 - 175/16384 * E2^4 -       411/65536 * E2^5);
F_2  =          3/8 * E2 + 3/16 * E2^2 + 213/2048 * E2^3 +  255/4096 * E2^4 +   166479/655360 * E2^5;
F_4  =                   21/256 * E2^2 +   21/256 * E2^3 +  533/8192 * E2^4 -   120563/327680 * E2^5;
% F_6  =                                   151/6144 * E2^3 +  155/4096 * E2^4 + 2767911/9175040 * E2^5;
F_6  =                                   151/6144 * E2^3 +  147/4096 * E2^4 + 2732071/9175040 * E2^5;
F_8  =                                                   1097/131072 * E2^4 -  273697/4587520 * E2^5;

Ge = G / e_0;
B   = Ge + F_2 * sin(2 * Ge) + F_4 * sin(4 * Ge) + F_6 * sin(6 * Ge) + F_8 * sin(8 * Ge);
