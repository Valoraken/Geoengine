function [X,Y] = Runge_Kutta3(fxy,x0,h,xmax,y0)
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
ti=x0;
yi=y0;
X=zeros((xmax-x0)/h,1);
Y=zeros((xmax-x0)/h,1);
% caculate Y for t=ti
for ti=x0:h:xmax
    X(i)=ti;
    Y(i)=yi;
    k1=fxy(ti,yi)
    k2=fxy(ti+0.5*h,yi+0.5*h*k1)
    k3=fxy(ti+h,yi-h*k1+2*h*k2)
    yi=yi+h/6*(k1+4*k2+k3)
    i=i+1;
end
