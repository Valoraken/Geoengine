% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Con_coe = Conformal_coe(B,m0,N,t,eta)
% This function is to calculate the coefficients of bivariate polynomials
% of degree 6.
% Input:  B   [rad]:latitude of the local origin P0(L0,B0) 
%         m0  [-] : scale factor
%         N   [m] : N = A/(1-E^2sin^2(B))^(1/2)
%         t   [-] : t = tan(B)
%         eta [-] : eta = (E^2/(1-E^2))cos(B)
%                 
% 
% Output: 
%         Con_coe: a 5*7 matrix contains all coefficients of bivariate polynomials
%                        
%   The index of all coefficients are 2 numbers [m n]
%   where m is refered to the exponent of b and n is refered to the exponent of l.

Con_coe = zeros(5,7);
% coefficients in polynomial of x
Con_coe(2,1)=m0*N*(1-eta^2+eta^4-eta^6);
Con_coe(3,1)=(3*m0/2)*N*t*(eta^2-2*eta^4);
Con_coe(4,1)=(m0/2)*N*(eta^2*(1-t^2)+eta^4*(-2+7*t^2));
Con_coe(5,1)=(-m0/2)*N*eta^2*t;

Con_coe(1,3)= m0/2*N*t*cos(B)*cos(B);
Con_coe(2,3)= m0/2*N*(1-t^2+eta^2*t^2-eta^4*t^2)*cos(B)*cos(B);
Con_coe(3,3)= m0/4*N*t*(-4+3*eta^2*(1-t^2))*cos(B)*cos(B);
Con_coe(4,3)= m0/3*N*(-1+t^2)*cos(B)*cos(B);
Con_coe(5,3)= m0/3*N*t*cos(B)*cos(B);

Con_coe(1,5)= m0/24*N*t*(5-t^2+9*eta^2)*cos(B)*cos(B)*cos(B)*cos(B);
Con_coe(2,5)= m0/24*N*(5-18*t^2+t^4+eta^2*(9-40*t^2-t^4))*(cos(B))^4;
Con_coe(3,5)= m0/6*N*t*(-7+5*t^2)*(cos(B))^4;

Con_coe(1,7)= m0/720*N*t*(61-58*t^2+t^4)*(cos(B))^6;
Con_coe(2,7)= m0/720*N*(61-479*t^2+179*t^4-t^6)*(cos(B))^6;

% coefficients in polynomial of y
Con_coe(1,2)= m0*N*cos(B);
Con_coe(2,2)= m0*N*t*(-1+eta^2-eta^4)*cos(B);
Con_coe(3,2)= m0/2*N*(-1+eta^2*(1-3*t^2)+eta^4*(-1+6*t^2))*cos(B);
Con_coe(4,2)= m0/6*N*t*(1+eta^2*(-10+3*t^2))*cos(B);
Con_coe(5,2)= m0/24*N*cos(B);

Con_coe(1,4)= m0/6*N*(1-t^2+eta^2)*(cos(B))^3;
Con_coe(2,4)= m0/6*N*t*(-5+t^2-eta^2*(4+t^2))*(cos(B))^3;
Con_coe(3,4)= m0/12*N*(-5+13*t^2+eta^2*(-4+8*t^2+3*t^4))*(cos(B))^3;
Con_coe(4,4)= m0/36*N*t*(41-13*t^2)*(cos(B))^3;

Con_coe(1,6)= m0/120*N*(5-18*t^2+t^4+eta^2*(14-58*t^2))*(cos(B))^5;
Con_coe(2,6)= m0/120*N*t*(-61+58*t^2-t^4)*(cos(B))^5;
Con_coe(3,6)= m0/240*N*(-61+418*t^2-121*t^4)*(cos(B))^5;
end

