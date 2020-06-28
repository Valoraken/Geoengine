% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Trb_clamda = Trb_c_lamda(B,m0,N,t,eta)
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
%         Trb_clamda: a 5*6 matrix contains all coefficients of bivariate polynomials
%                        
%   The index of all coefficients are 2 numbers [m n]
%   where m is refered to the exponent of x and n is refered to the exponent of y.
Trb_clamda=zeros(5,6);
% coefficients that compute c
Trb_clamda(1,2)= 1*t/(m0*N);
Trb_clamda(2,2)= (1)*(1+t^2+eta^2)/(m0^2*N^2);
Trb_clamda(3,2)= (1*t)*(1+t^2-eta^2)/(m0^3*N^3);
Trb_clamda(4,2)= (1)*(1+4*t^2+3*t^4)/(3*(m0^4*N^4));
Trb_clamda(5,2)= (1*t)*(2+5*t^2+3*t^4)/(3*(m0^5*N^5));

Trb_clamda(1,4)= -(1*t)*(1+t^2-eta^2)/(3*(m0^3*N^3));
Trb_clamda(2,4)= -(1)*(1+4*t^2+3*t^4)/(3*(m0^4*N^4));
Trb_clamda(3,4)= -(2*t)*(2+5*t^2+3*t^4)/(3*(m0^5*N^5));
Trb_clamda(4,4)= -(2)*(2+17*t^2+30*t^4+15*t^6)/(9*(m0^6*N^6));

Trb_clamda(1,6)= (1*t)*(2+5*t^2+3*t^4)/(15*(m0^5*N^5));
Trb_clamda(2,6)= (1)*(2+17*t^2+30*t^4+15*t^6)/(15*(m0^6*N^6));
Trb_clamda(3,6)= (1*t)*(17+77*t^2+105*t^4+45*t^6)/(15*(m0^7*N^7));

% coefficients that compute lamda
Trb_clamda(1,3) = (1)*(1+eta^2)/(2*m0*N^2);
Trb_clamda(2,3) = -2*t*eta^2/(m0^2*N^3);
Trb_clamda(3,3) = 1*eta^2*(-1+t^2)/(m0^3*N^4);

Trb_clamda(1,5) = 1*(1+6*eta^2)/(24*m0^3*N^4);
end