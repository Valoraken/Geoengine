clear all
close all
clc
load AnscombeQuartet.mat
xhat_1=zeros(2,4);
xhat_2=zeros(3,4);
yhat_1=zeros(11,4);
yhat_2=zeros(11,4);
r1=11-2;
r2=11-3;
for i=1:4
    i
A=ones(11,1);
A(:,2)=data(:,2*i-1)
N=A'*A
inv(N)
xhat=inv(N)*A'*data(:,2*i)
yhat=A*xhat
ehat=data(:,2*i)-yhat
xhat_1(:,i)=xhat;
yhat_1(:,i)=yhat;
ehat_1(:,i)=ehat;
sumehat_1(i)=ehat'*ehat
nsumehat_1(i)=ehat'*ehat/r1
p=0:20;
y_1=xhat_1(2,i)*p+xhat_1(1,i);
subplot(2,4,2*i-1)
scatter(data(:,2*i-1),data(:,2*i));hold on
plot(p,y_1)
xlabel('x'); ylabel('y');
title({['Data set ',num2str(i)],['Simple linear regression']})
grid on

A(:,3)=data(:,2*i-1).^2
N=A'*A
inv(N)
xhat=inv(N)*A'*data(:,2*i)
yhat=A*xhat
ehat=data(:,2*i)-yhat
xhat_2(:,i)=xhat;
yhat_2(:,i)=yhat;
ehat_2(:,i)=ehat;
sumehat_2(i)=ehat'*ehat
nsumehat_2(i)=ehat'*ehat/r2
y_2=xhat_2(1,i)+xhat_2(2,i)*p+xhat_2(3,i)*p.^2;
subplot(2,4,2*i)
scatter(data(:,2*i-1),data(:,2*i));hold on
plot(p,y_2)
xlabel('x'); ylabel('y');
title({['Data set ',num2str(i)],['polynomial regression']})
grid on
end
print('-depsc2', 'figure.eps');
print('-dpng', '-r800', 'figure.png')

