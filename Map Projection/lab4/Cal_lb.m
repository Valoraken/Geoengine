% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Ell_points = Cal_lb(m0,x,y,B)
% This function is to compute the ellipsoidal corrdinates,meridian convergance and distortion of a given
% points(x,y) in the Cartesian coordinates
% Input: 
%       m0  [-] : scale factor
%       x   [m] :x coordinate of the given point(x,y) 
%       y   [m] :y coordinate of the given point(x,y) 
%       B   [rad]:latitude of the local origin P0(L0,B0) 
% Output:
%       Ell_points = [b,l,c,lamda]
%       b   [rad]: latitude  differences between the given point P(Lp,Bp) and the local origin P0(L0,B0)
%       l   [rad]: longtitude differences between the given point P(Lp,Bp) and the local origin P0(L0,B0)  
%       c   [rad]: maridian convergance
%       lamda [-]: distortion


A_B= 6377397.155;
E_2B=0.0066743722;

%% compute M N t and eta at the point B0=B
N = A_B/(1-E_2B*sin(B)^2)^0.5;
t = tan(B);
eta = cos(B)*(E_2B/(1-E_2B))^0.5;

%% compute the coefficients of x and y also the meridian convergence and distortion

Trb_coe = Transback_coe(B,m0,N,t,eta);
Trb_clamda = Trb_c_lamda(B,m0,N,t,eta);
b= [1 x x^2 x^3 x^4]*Trb_coe(:,1)+y^2.*[1 x x^2 x^3 x^4]*Trb_coe(:,3)+y^4.*[1 x x^2 x^3 x^4]*Trb_coe(:,5)+y^6.*[1 x x^2 x^3 x^4]*Trb_coe(:,7);
l = y.*[1 x x^2 x^3 x^4]*Trb_coe(:,2)+ y^3.*[1 x x^2 x^3 x^4]*Trb_coe(:,4)+y^5.*[1 x x^2 x^3 x^4]*Trb_coe(:,6)+y^7.*[1 x x^2 x^3 x^4]*Trb_coe(:,8);
c = y.*[1 x x^2 x^3 x^4]*Trb_clamda(:,2)+ y^3.*[1 x x^2 x^3 x^4]*Trb_clamda(:,4)+ y^5.*[1 x x^2 x^3 x^4]*Trb_clamda(:,6);
lamda = m0+ y^2.*[1 x x^2 x^3 x^4]*Trb_clamda(:,3)+ y^4.*[1 x x^2 x^3 x^4]*Trb_clamda(:,5);

Ell_points=[b,l,c,lamda];

end