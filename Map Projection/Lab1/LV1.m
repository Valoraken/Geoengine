%Xiao, Tianqi 3371477
%lab1 02.11.18
% Plot the surface:
%    X(U,V)=A/Bsinh(V)cos(U)*E1*+A/Bsinh(V)sin(U)*E2*+AU*E3*
% -2*pi<=U<=2*pi, 0<=V<=6
% A ,B=const, Assume A=B=3
%
syms X Y Z U V;
X= @(U,V) sinh(V).*cos(U);
Y= @(U,V) sinh(V).*sin(U);
Z= @(U,V) 3*U;
%Plot the parameterized mesh
fmesh(X,Y,Z,[-2*pi 2*pi 0 6])
xlabel('X');ylabel('Y');zlabel('Z')
title ({'"Helicoid 2"';' [X= sinh(V)*cos(U), Y= sinh(V)*sin(U),  Z=3*U] '})
alpha(0.7)
