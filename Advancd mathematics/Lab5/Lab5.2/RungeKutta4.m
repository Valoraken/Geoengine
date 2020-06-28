% Xiao, Tianqi 3371477
% Lab 5 27.11.18

function [X,Y] = RungeKutta4(fxy,x0,h,xmax,y0)
% This function is to implement the Runge-Kutta method of order four for numerical integration

% INPUT:
%  fxy is a function handle which provide the differential equation
%  x0 is the initial value of x
%  h is the stepwidth
%  xmax is the maximum value of x 
%  y0 is the initial value of y

% OUTPUT:
% Y is the value of the solution function at X
% X is the value of x in the interval of [0,6] with the stepwidth h
% the result will also be visualized

%%   check input
if class(fxy)~='function_handle'
    error('incorrect fucntion type! ');
end

if length(x0)~=1 
    error('incorrect dimension of x0');
end
if length(h)~=1
     error('incorrect dimension of h');
end
if length(xmax)~=1
     error('incorrect dimension of xmax');
end

if iscolumn(y0)==0
    error('y0 must be scalar or column vector ')
end

if size(y0)~=size(fxy)
    error('The dimention of fxy and y0 must be the same')
end

%% implement the Runge-Kutta method and visualize the result
%inicialize
i=1;
xi=x0;
yi=y0;
X=zeros((xmax-x0)/h,1);
Y=zeros((xmax-x0)/h,1);
% caculate Y for X=xi
for xi=x0:h:xmax
    X(i)=xi;
    Y(i)=yi(1);
    k1=fxy(xi,yi);
    k2=fxy(xi+0.5*h,yi+0.5*h*k1);
    k3=fxy(xi+0.5*h,yi+0.5*h*k2);
    k4=fxy(xi+h,yi+h*k3);
    yi=yi+h/6*(k1+2*k2+2*k3+k4);
    i=i+1;
end
% visualize
plot(X,Y)
