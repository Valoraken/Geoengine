% Physical Geodesy Exercise 2
% NAME: Xiao, Tianqi, Wang Zhenqiao
% Student ID: 3371477, 3371590
% 23.5.2019

function [a] = a_shell(R_in,R_out,rho,r)
% This function is to calculate the attraction of a
% spherical shell depending on the inner and outer radius
%    INPUT  :
%   R_in  : inner radius of the spherical shell
%   R_out : outer radius of the spherical shell
%    rho  : density of the solid sphere
%    r    : radial coordinate r of the point P.
%    OUTPUT :
%    a   : attraction at the point P
G= 6.672e-11;
% inner part
a(r<=R_in) = 0;
% in shell
a(r>R_in & r<=R_out) = -4/3*pi*G*rho*(r(r>R_in & r<=R_out).^3-R_in^3)./r(r>R_in & r<=R_out).^2;
% outer part
a(r>R_out) = -4/3*pi*G*rho*(R_out^3-R_in^3)./r(r>R_out).^2;
end


