% Lab 1 Task 2 Relative gravimetry
% Group 2
% 12.05.2019
clear all
clc

%% Measured data and Pre-Caculation
% Given values
delta_h =[6.897 17.260]';                                                           % Height differences dh(13) dh(16) 
sigma_dh =0.03;                                                                       % Presition of Height Differences 
GM =3.986004415e14;                                                                 % Geocentric gravity constant
R =6378136                                                                          % Radius of the Earth

% Caculated data from the Excel
dy =[-1.8905 1.893 -1.884 -2.962 2.9695 -2.971]';                                   % Differences between Observations [mGal]
Qdy =[0.1550445 0.17860425 0.058785 0.03541125 0.2170535 0.2165355]';               % Covariance matrix of dy [mGal^2]
Qy =[0.01374025 0.14130425 0.0373 0.021485 0.01392625 0.20312725 0.01340825]';      % Covariance matrix of y  [mGal^2]
dt =[432 403 309 436 307 296]';                                                     % Average time difference of each station [s]

% Define Design matrix
fprintf('Design Matrix\r\n');
D =[1  0; -1  0;1  0;0  1;0 -1;0  1];
A =[D dt]
% unknowns: delta_g1_3 delta_g3_6 drift
% Define Weight matrix
fprintf('Weight Matrix\r\n');
P = diag(1./Qdy)
%% Least-Square Adjustment
N =A'*P*A;
x_hat =N\A'*P*dy;
dy_hat =A*x_hat;                                          
e_hat =dy-dy_hat;

% Orthogonality check
check =A'*P*e_hat;
if norm(check) > 1e-10
   error('Orthogonality check failure...\r\n');
else 
    fprintf('Orthogonality check pass...\r\n');
end

fprintf('Solution£º\r\n')
fprintf('  dg(13) = %.4f mGal\r\n  dg(36) = %.4f mGal\r\n  dg(16) = %.4f mGal\r\n',x_hat(1),x_hat(2),x_hat(1)+x_hat(2))
%% Accuracy
L = inv(N)*A'*P;
Q = diag(Qdy); 
sigma2_x_hat = L*Q*L';

for i=1:6
sig_dy(i)=sqrt(Qy(i)+Qy(i+1));
end

%% Gradient
delta_g=[x_hat(1) x_hat(1)+x_hat(2)]';

grad = -delta_g./delta_h;
fprintf('Gravity Verticle Gradient:\r\n  grad_13= %.4f mGal/m\r\n  grad_16= %.4f mGal/m\r\n',grad)

% accuracy
sigma2_dg=[sigma2_x_hat(1,1) sigma2_x_hat(1,1)+sigma2_x_hat(2,2)]';
sigma_grad2=(1./(delta_h)).^2.*sigma2_dg+(delta_g./delta_h.^2).^2.*sigma_dh^2;
fprintf('Accuracy of Gravity Verticle Gradient:\r\n  grad_13= %.4s\r\n  grad_16= %.4s\r\n',sigma_grad2)
% Calculated gradient
V=GM/R;
Grad_cal=2*V/R^2*1e5;
fprintf('Calculated Gravity Verticle Gradient:\r\n  %.4f mGal/m\n',Grad_cal)