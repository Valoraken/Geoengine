% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
function Pro_points = Cal_xy(m0,l,b,B)
% This function is to compute the Cartisian corrdinates,meridian convergance and distortion of a given
% points(Lp,Bp) in the GRS80 ellipsoidal coordinates
% Input: 
%       m0  [-] : scale factor
%       l   [rad]:longtitude differences between the given point P(Lp,Bp) and the local origin P0(L0,B0)  
%       b   [rad]:latitude  differences between the given point P(Lp,Bp) and the local origin P0(L0,B0)
%       B   [rad]:latitude of the local origin P0(L0,B0) 
% Output:
%       Pro_points = [x,y,c,lamda]
%       (x,y) [m]: the cartiesian coordinates of the given point in the map
%       c   [rad]: maridian convergance
%       lamda [-]: distortion


% GRS80 ellipsoid
A= 6378137;
E_2 =0.00669438002;
%% compute M N t and eta at the point B0=B
N = A/(1-E_2*sin(B)^2)^0.5;
t = tan(B);
eta = cos(B)*(E_2/(1-E_2))^0.5;

%% compute the coefficients of x and y also the meridian convergence and distortion

Con_coe = Conformal_coe(B,m0,N,t,eta);
Coe_clamda = Coe_c_lamda(B,m0,N,t,eta);
x= [1 b b^2 b^3 b^4]*Con_coe(:,1)+l^2.*[1 b b^2 b^3 b^4]*Con_coe(:,3)+l^4.*[1 b b^2 b^3 b^4]*Con_coe(:,5)+l^6.*[1 b b^2 b^3 b^4]*Con_coe(:,7);
y = l.*[1 b b^2 b^3 b^4]*Con_coe(:,2)+ l^3.*[1 b b^2 b^3 b^4]*Con_coe(:,4)+l^5.*[1 b b^2 b^3 b^4]*Con_coe(:,6);
c = l.*[1 b b^2 b^3]*Coe_clamda(:,2)+ l^3.*[1 b b^2 b^3]*Coe_clamda(:,4)+ l^5.*[1 b b^2 b^3]*Coe_clamda(:,6);
lamda =m0+ l^2.*[1 b b^2 b^3]*Coe_clamda(:,3)+ l^4.*[1 b b^2 b^3]*Coe_clamda(:,5)+ l^6.*[1 b b^2 b^3]*Coe_clamda(:,7);

Pro_points=[x,y,c,lamda];
end

