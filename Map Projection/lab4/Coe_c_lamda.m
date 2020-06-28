% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Coe_clamda = Coe_c_lamda(B,m0,N,t,eta)
 
% This function is to calculate the coefficients of the polinomialsof degree 6 that
% compute the meridian convergence c and the distortion("scale factor")
% 
% Input:  B   [rad]:latitude of the local origin P0(L0,B0) 
%         m0  [-] : scale factor
%         N   [m] : N = A/(1-E^2sin^2(B))^(1/2)
%         t   [-] : t = tan(B)
%         eta [-] : eta = (E^2/(1-E^2))cos(B)
%                 
% 
% Output: 
%         Coe_clamda: a 4*7 matrix contains all coefficients of bivariate polynomials
%                        
%   The index of all coefficients are 2 numbers [m n]
%   where m is refered to the exponent of b and n is refered to the exponent of l.

Coe_clamda = zeros(4,7);

% coefficients of the polynomials that compute c
Coe_clamda(1,2)= t*cos(B);
Coe_clamda(2,2)= cos(B);
Coe_clamda(3,2)= -t/2*cos(B);
Coe_clamda(4,2)= -1/6*cos(B);

Coe_clamda(1,4)= 1/3*t*(1+3*eta^2)*(cos(B))^3;
Coe_clamda(2,4)= 1/3*(1-2*t^2+eta^2*(3-12*t^2))*(cos(B))^3;
Coe_clamda(3,4)= 1/6*t*(-7+2*t^2)*(cos(B))^3;

Coe_clamda(1,6)= 1/15*t*(2-t^2)*(cos(B))^5;

% coefficients of the polynomials that compute lamda

Coe_clamda(1,3)= m0/2*(1+eta^2)*(cos(B))^2;
Coe_clamda(2,3)= -m0*t*(1+2*eta^2)*(cos(B))^2;
Coe_clamda(3,3)= m0/2*(-1+t^2+eta^2*(-2+6*t^2))*(cos(B))^2;
Coe_clamda(4,3)= 2*m0/3*t*(cos(B))^2;

Coe_clamda(1,5)= m0/24*(5-4*t^2+eta^2*(14-28*t^2))*(cos(B))^4;
Coe_clamda(2,5) = m0/6*t*(-7+2*t^2)*(cos(B))^4;

Coe_clamda(1,7)= m0/720*(61-148*t^2+16*t^4)*(cos(B))^6;
end

