% Xiao,Tianqi
% Lab 1 Task 2
% 21.04.2019
% the ODE:  y'=-2ty/(y^2-t^2)
% interval of t: [0,4]
% stepwidth: h=1
% initial values: t0=0 y0=4
%
clear 
clc
% the rewritten function handle
f=@(t,y) -2*t*y/(y^2-t^2);
t0=0;
h=0.01;
tmax=4;
y0=4;
% visualize
[X,Y]= RungeKutta3(f,t0,h,tmax,y0);
plot(X,Y,'r')
fprintf('y(4) by 3rd order method: %f \r',Y(tmax/h+1))
hold on
[x,y]=RungeKutta4(f,t0,h,tmax,y0);
plot(x,y,'b')
grid on
fprintf('y(4) by 4th order method: %f \r',y(tmax/h+1))
legend('3rd order','4th order')
xlabel('t');
ylabel('y');
title({'Solution of y¡¯=-2ty/(y^2-t^2) ';('t0=0  y0=4')})
