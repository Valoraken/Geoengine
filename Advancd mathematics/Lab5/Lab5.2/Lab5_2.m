% Xiao, Tianqi 3371477
% Lab 5-2 27.11.18

% the ODE:  y''+xy=0
% interval of x: [0,6]
% stepwidth: 0.01
% initial values: x0=0 y0= 0.35502805388
%
% the ODE is second order, rewrite into the system:
% Y'=[0  1] Y ,   with Y=[y,y']'
%    [-x 0] 
clear 
clc
% the rewritten function handle
fxy=@(x,y)[y(2);-x*y(1)];

x0=0; h=0.01; xmax=6;
% no y0' given assume y0'=1
y0=[0.35502805388;1];

[X,Y]= RungeKutta4(fxy,x0,h,xmax,y0);
title({'Solution of y"+ xy = 0 ';('x0=0  y0=0.3502805388  y0''=1')})