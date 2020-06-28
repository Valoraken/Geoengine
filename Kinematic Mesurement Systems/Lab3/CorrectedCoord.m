function [X_corr, Y_corr] = CorrectedCoord(T, Hz, V, S, detat)
% [X_corr, Y_corr] = CorrectedCoord(T, Hz, V, S, detat)
% 
% CORRECTEDCOORD aims to compute the corrected coordinates of tracking targets
% 
% n = length(T);
% Input:
% T:       Time stamp              (column vector)        [n*1]
% Hz:      Horizontal angle        (column vector)        [n*1]
% V:       Vertical angle          (column vector)        [n*1]
% S:       Slope distance          (column vector)        [n*1]
% detat:   sychronization error    (scalar)
% Output:
% X_corr: Corrected X coordinates  (column vector)        [n*1]
% Y_corr: Corrected Y coordinates  (column vector)        [n*1]
% Name: YI HONG
% Matriculation number:3294211
%
n = length(Hz);
% initialize vectors
X_corr = zeros(n, 1);                                         
Y_corr = zeros(n, 1);

% convert slope distance to horizontal distance
V = deg2rad(V*0.9);      
S_H = S.*sin(V);    

% real distance
dt = zeros(n, 1);
detas = zeros(n, 1);
for i = 1: n-1
    dt(i+1) = (T(i+1)-T(i))/1000;
    detas(i+1) = detat*(S_H(i)-S_H(i+1))./dt(i+1);
end
S_real = S_H + detas;

% Corrected coordinates
Hz = deg2rad(Hz*0.9);
Xs = 4.71174;                                                % position of the tachymeter
Ys = 2.2536;
Orient = 2.10434;                                            % orientation angle [rad]
X_corr = Xs + S_real.*cos(Orient+Hz);
Y_corr = Ys + S_real.*sin(Orient+Hz);
end