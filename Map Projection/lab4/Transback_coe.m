% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Trb_coe = Transback_coe(B,m0,N,t,eta)
% This function is to calculate the coefficients of transback polynomials of degree 6.
% Input:  B   [rad]:latitude of the local origin P0(L0,B0) 
%         m0  [-] : scale factor
%         N   [m] : N = A/(1-E^2sin^2(B))^(1/2)
%         t   [-] : t = tan(B)
%         eta [-] : eta = (E^2/(1-E^2))cos(B)
%                 
% 
% Output: 
%         Trb_coe: a 5*8 matrix contains all coefficients of transback polynomials
%                        
%   The index of all coefficients are 2 numbers [m n]
%   where m is refered to the exponent of x and n is refered to the exponent of y.
Trb_coe=zeros(5,8);
% coefficients of b(latitude)
Trb_coe(2,1)= (1)*(1+eta^2)/(m0*N);
Trb_coe(3,1)= -(3*t)*(eta^2+eta^4)/(2*(m0^2*N^2));
Trb_coe(4,1)= (1)*(eta^2*(-1+t^2)+eta^4*(-2+6*t^2))/(2*(m0^3*N^3));
Trb_coe(5,1)= (1)*t*eta^2/(2*(m0^4*N^4));

Trb_coe(1,3)= -(1*t)*(1+eta^2)/(2*(m0^2*N^2));
Trb_coe(2,3)= (1)*(-1-t^2+eta^2*(-2+2*t^2)+eta^4*(-1+3*t^2))/(2*(m0^3*N^3));
Trb_coe(3,3)= (1*t)*(-2-2*t^2+eta^2*(9+t^2))/(4*(m0^4*N^4));
Trb_coe(4,3)= (1)*(-2-8*t^2-6*t^4+eta^2*(7-6*t^2+3*t^4))/(12*(m0^5*N^5));
Trb_coe(5,3)= -(1*t)*(2+5*t^2+3*t^4)/(6*(m0^6*N^6));

Trb_coe(1,5)= (1*t)*(5+3*t^2+eta^2*(6-6*t^2))/(24*(m0^4*N^4));
Trb_coe(2,5)= (1)*(5+14*t^2+9*t^4+eta^2*(11-30*t^2-9*t^4))/(24*(m0^5*N^5));
Trb_coe(3,5)= (1*t)*(7+16*t^2+9*t^4)/(12*(m0^6*N^6));
Trb_coe(4,5)= (1)*(7+55*t^2+93*t^4+45*t^6)/(36*(m0^7*N^7));

Trb_coe(1,7)= -(1*t)*(61+90*t^2+45*t^4)/(720*(m0^6*N^6));
Trb_coe(2,7)= -(1)*(61+331*t^2+495*t^4+225*t^6)/(720*(m0^7*N^7));

% coeffcients of l (longtitude)
Trb_coe(1,2)= 1/(m0*N*cos(B));
Trb_coe(2,2)= 1*t/(m0^2*N^2*cos(B));
Trb_coe(3,2)= (1)*(1+2*t^2+eta^2)/(2*(m0^3*N^3*cos(B)));
Trb_coe(4,2)= (1*t)*(5+6*t^2+eta^2)/(6*(m0^4*N^4*cos(B)));
Trb_coe(5,2)= (1)*(5+28*t^2+24*t^4)/(24*(m0^5*N^5*cos(B)));

Trb_coe(1,4)= -(1)*(1+2*t^2+eta^2)/(6*(m0^3*N^3*cos(B)));
Trb_coe(2,4)= -(1*t)*(5+6*t^2+eta^2)/(6*(m0^4*N^4*cos(B)));
Trb_coe(3,4)= -(1)*(5+28*t^2+24*t^4+eta^2*(6+8*t^2))/(12*(m0^5*N^5*cos(B)));
Trb_coe(4,4)= -(1*t)*(61+180*t^2+120*t^4)/(36*(m0^6*N^6*cos(B)));

Trb_coe(1,6)= (1)*(5+28*t^2+24*t^4+eta^2*(6+8*t^2))/(120*(m0^5*N^5*cos(B)));
Trb_coe(2,6)= (1*t)*(61+180*t^2+120*t^4)/(120*(m0^6*N^6*cos(B)));
Trb_coe(3,6)= (1)*(61+662*t^2+1320*t^4+720*t^6)/(240*(m0^7*N^7*cos(B)));

Trb_coe(1,8)= -(1)*(61+662*t^2+1320*t^4+720*t^6)/(5040*(m0^7*N^7*cos(B)));

end

