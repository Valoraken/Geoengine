% Lab 1 Task 1 Absolute gravimetry
% Group 2
% 12.05.2019
clear all
clc
%% measured data
T= 113.79;                          % Total time
n=15;                               % Number of periods
l0=0.202;                            % Distance between the mass centre and the ground
L=14.347;                           % Distance between the privot point and the 1st floor

%% caculation of absolute gravity
t=T/n;                              % Period
l=L-l0;                             % Length of the pendulum
g=l*(2*pi/t)^2                      % Gravity

%% accuracy of gravity
sigma_t=0.001;                      % Accuracy of time
sigma_l=0.001;                      % Accuracy of length

% Error propagation
sigma_g2=((2*pi/t)^2)^2*(sigma_l)^2+((-2*l/t)*((2*pi/t)^2))^2*(sigma_t)^2;
sigma_g=sqrt(sigma_g2)              % Accuracy of gravity
%% accuracy of time measurement required
sigma_g0= 1.0e-5;                   % Accuracy of gravity required

sigma_t0=(sigma_g0-(2*pi/t)^2*sigma_l)/(-2*l/t)*((2*pi/t)^2)       % Accuracy of time required
