close all
clear 
clc

h =21.86;                           % height of the tower [m]
x0 =13                              % horizontal distance [m]
% observations
t = [0 3.2 7.5 12.1 17.8];            % time [s]

z = [36.25 55.06 76.66 104.42];     % measurement distance[m]
R =0.01;                            % mesurement noise [m^2]
% z = (x^2+h^2)^0.5
H = x0/sqrt(x0^2+h^2);

% random walk: x'=0+w(t)
F = 0;
sigma_2 =4;                               % varience of process noise
sigma0 =1; 
%% Apply EKF
Pnn = sigma0;
xnn = x0;
x_pre=[];
x_fwd=[];
P_fwd=[];
I = eye(size(sigma0));
for n=1:4
    dt=t(n+1)-t(n);
    phi=exp(F*dt);
    Q = sigma_2*dt;
    xnn_p= phi*xnn;
    Pnn_p= phi*Pnn*phi'+Q
    H =  xnn_p/sqrt(xnn_p^2+h^2);
    K= Pnn_p*H'*(H*Pnn_p*H'+ R)^(-1);
    xnn=xnn_p+K*(z(n)-sqrt(xnn_p^2+h^2));
    Pnn=(I-K*H)*Pnn_p;
    
    % save result 
    x_pre=[x_pre xnn_p];
    x_fwd=[x_fwd xnn];
    Sigma_fwd(:,n)=sqrt(diag(Pnn));
    P_fwd=[P_fwd Pnn];
end
figure
hold on

plot(t(2:5),x_fwd,'*-');
plot(t(2:5),sqrt(z.^2-h^2),'+-');
legend('x EKF','x Observe')
xlabel('t [s]');
ylabel('x [m]');
figure
hold on
plot(t(2:5),Sigma_fwd);
xlabel('t [s]');
ylabel('\sigmax(t) [m]');