 % Lab 4 part 2
% Xiao, Tianqi 3371477
% 08.01.19
clear
x=[-4 -3 -2 1 2 4]';
y=[-26 -775 -1104 -1011 -980 702]';
%% degree of 3
% numbers of observations n = 6
% numbers of unknowns k = 3 + 1 = 4 
n=6;
k=4;
r=n-k;
% initialize design matrix A
A_1=zeros(n,k);
% build A
for i =1:k
A_1(:,i)=x.^(i-1);
end
% A_1
% compute coefficents
xhat_1=A_1\y
yhat_1=A_1*xhat_1
ehat_1=y-yhat_1
sumehat_1=ehat_1'*ehat_1
nsumehat_1=ehat_1'*ehat_1/r
% plot the polynomial
p=-5:0.1:5;
Y_1=xhat_1(1)+xhat_1(2).*p+xhat_1(3).*p.^2+xhat_1(4).*p.^3;
subplot(3,1,1)
scatter(x,y);hold on
plot(p,Y_1)
xlabel('x'); ylabel('y');
text(0,500,['eTe= ',num2str(sumehat_1)])
title({['y=',num2str(xhat_1(1)),'+',num2str(xhat_1(2)),'x+',num2str(xhat_1(3)),'x^2+',num2str(xhat_1(4)),'x^3'],'Polynomial of degree 3'})
grid on

%% degree of 4
% numbers of observations n = 6
% numbers of unknowns k = 4 + 1 = 5 
n=6;
k=5;
r=n-k;
% initialize design matrix A
A_2=zeros(n,k);
% build A
for i =1:k
A_2(:,i)=x.^(i-1);
end
% A_2
% compute coefficents
xhat_2=A_2\y
yhat_2=A_2*xhat_2
ehat_2=y-yhat_2
sumehat_2=ehat_2'*ehat_2;
nsumehat_2=ehat_2'*ehat_2/r;
p=-5:0.1:5;
Y_2=xhat_2(1)+xhat_2(2).*p+xhat_2(3).*p.^2+xhat_2(4).*p.^3+xhat_2(5).*p.^4;
subplot(3,1,2)
scatter(x,y);hold on
plot(p,Y_2)
xlabel('x'); ylabel('y');
text(0,500,['eTe= ',num2str(sumehat_2)])
title({['y=',num2str(xhat_2(1)),'+',num2str(xhat_2(2)),'x+',num2str(xhat_2(3)),'x^2+',num2str(xhat_2(4)),'x^3+',num2str(xhat_2(5)),'x^4'],['Polynomial of degree 4']})
grid on

%% degree of 5 
% numbers of observations n = 6
% numbers of unknowns k = 5 + 1 = 6 
n=6;
k=6;
r=n-k;
% initialize design matrix A
A_b=zeros(n,k);
% build A
for i =1:k
A_b(:,i)=x.^(i-1);
end
% A_b
% compute coefficents
xhat_b=A_b\y
yhat_b=A_b*xhat_b
ehat_b=y-yhat_b
sumehat_b=ehat_b'*ehat_b;
nsumehat_b=ehat_b'*ehat_b/r;

%plot the polynomial
p=-5:0.1:5;
Y_b=xhat_b(1)+xhat_b(2).*p+xhat_b(3).*p.^2+xhat_b(4).*p.^3+xhat_b(5).*p.^4+xhat_b(6).*p.^5;
subplot(3,1,3)
scatter(x,y);hold on
plot(p,Y_b)
xlabel('x'); ylabel('y');
text(0,500,['eTe= ',num2str(sumehat_b)])
title({['y=',num2str(xhat_b(1)),'+',num2str(xhat_b(2)),'x+',num2str(xhat_b(3)),'x^2+',num2str(xhat_b(4)),'x^3+',num2str(xhat_b(5)),'x^4+',num2str(xhat_b(6)),'x^5'],'Polynomial of degree 5'})
fprintf (['y=',num2str(xhat_b(1)),'+',num2str(xhat_b(2)),'x+',num2str(xhat_b(3)),'x^2+',num2str(xhat_b(4)),'x^3+',num2str(xhat_b(5)),'x^4+',num2str(xhat_b(6)),'x^5','Polynomial of degree 5'])
grid on


%% compute Li(x) 

polyLagrange(x,y)
