clear all
close all
clc
omega0 = 0.03;
f = 0.99999;
deltaT = 1;


F = [0,1;-omega0.^2,-2*f*omega0];
phi=expm(F*deltaT);
x = zeros(2,100001);

for i = 1:100000
    w = normrnd(0,1);
     x(:,i+1) = phi*x(:,i)+ [0;w];
end

tau = -1000:1000;
sigma=var(x(1,:));
beta=omega0*sqrt(1-f^2);
% define covarience
COVx1x1=sigma*exp(-f*omega0*abs(tau)).*(cos(beta*abs(tau))+(f*omega0/beta)*sin(beta*abs(tau)));
figure(1)
plot(tau,xcov(x(1,:),x(1,:),1000,'biased'),tau,COVx1x1);
xlabel('tau');ylabel('Covx1x1')
legend('empiric corarience ', 'theoretic covarience')
