% Xiao, Tianqi 3371477
% Lab 5-4 27.11.18
clear
clc
% test matrix 
A=[9,-3,2,1,2;-3,2,1,8,2;2,1,-4,-4,1;1,8,-4,-4,1;2,2,1,1,10];
% max iterations
iterMax=100;
% tolerance
tol=1e-10;
% result
A=jacobipq(A,iterMax,tol)