% Xiao Tianqi 3371477
% Lab 1 26.10.18
% The differential equation  2y''''+4y'''-34y''-36y'+144y = 0 is a Forth-order differential equation. 
% So its characteristic equation 2*r^4+4*r^3-34*r^2-36r+144=0 has four roots.

syms x C1 C2 C3 C4
% the coefficient vector an:
an=[144,-36,-34,4,2];
% all coefficient can be divided by 2
% So guess the solution of the characteristic equations is the divisor of 72
x0=[-72,-36,-24,-18,-12,-9,-8,-6,-4,-3,-2,-1,0,1,2,3,4,6,8,9,12,18,24,36,72];

% Implement the Horner scheme to find the roots
for i=1:length(x0)
    y(i)=horner(an,x0(i));
end

root=x0(find(y==0));
fprintf('The roots of characteristic equation are: ');
fprintf('%d ',root);
% Accroding to the root, the characteristic equation has four real roots:-4, -3, 2, 3
% Hence the general solution of the differential equation is
%       y = C1*exp^(-4x)+C2*exp^(-3x)+C3*exp^(2x)+C4*exp^(3x)
fprintf('\nThe general solution of the differential equation is: y= C1*exp^(%dx)+C2*exp^(%dx)+C3*exp^(%dx)+C4*exp^(%dx)\n',root );

