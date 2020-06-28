% Task 3
close all
clear;
clc;
dt=0.25;
F=[0,-0.3,0.7;0.3,0,0;-0.7,0,0];
Phi=expm(F*dt);
x=[2.3,1.7,6.4]';
solution=zeros(3,40);
solution(:,1)=x;
for i=1:39
    x=Phi*x;
    solution(:,i+1)=x;
end
t=0.25:.25:10;
plot(t,solution(1,:)); 
hold on;
plot(t,solution(2,:)); 
plot(t,solution(3,:));
xlabel('t');ylabel('x')
legend('x1','x2','x3');