% Physical Geodesy Lab 1
% Absolute and relative gravimetry
% Written by group 3
%
%% Data preparation
clc
clear all
dy = [-3.2880 3.2831 -1.8870 -2.9575 2.9710 -2.9695]';                       % Observation difference dy[mGal]
dt = [831 660 484 536 440 416]';                                             % Average time difference of each station dt[s]
Qdy = [0.01715504 0.03889565 0.09700050 0.08054550 0.04936625 0.05401625]';  % Accuracy square of Qdy[mGal^2]
                                                                             % Parameters above obtained from Excel
D = [ 1  0;                                                                   
     -1  0;  
      1  0;  
      0  1;  
      0 -1;  
      0  1];
A = [D dt];                                                                   % Design matrix A
P = diag(1./Qdy);                                                             % Weight matrix P

%% Least-square adjustment
N = A'*P*A;
x_hat = N\A'*P*dy
dy_hat = A * x_hat;                                          
e_hat = dy - dy_hat ;                                                                        
check = A'*P*e_hat;                                                          % Orthogonality check
if norm(check) > 1e-10
   error('Orthogonality check failure');
else 
    fprintf('Orthogonality check pass');
end

%% Error propogation
L = inv(N)*A'*P;
Q = diag(Qdy); 
sigma2_x_hat = L*Q*L';

