% Dynamic System Estimation
% 12/5/2019
% task 2
close all
clear all
clc

%% load data and initial values
x_out=zeros(20,100);
W=load('random02.txt');
dt=1;
F=[0,1;0,0];
Phi=expm(F*dt);
Phi
%% compute the stochastic process for 100 times with 20 steps each
N=20;
for j=1:100
w=W(:,j);
x=[0 0]';
x_solution=zeros(2,N);
x_solution(:,1)=x;

for  i=1:20
    x=Phi*x+[0;w(i)];
    x_solution(:,i)=x;
end
x_out(:,j)=x_solution(1,:);
end



%% compute the mean value and its variance
Mean=zeros(20:1);
Varience=zeros(20:1);

for i=1:20
Mean(i)=mean(x_out(i,:));
end
Varience = var(x_out,1,2);
plot(Mean);hold on
plot(Varience);
legend('mean','varience')
Mean'
Varience



