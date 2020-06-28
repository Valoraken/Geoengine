clear


%% Load data
data=load('KF_task2.txt');
% time
t=data(:,1);       % time
% observation
z=data(:,2);       % position x1

%% Model : m * x"(t) + b * x'(t) + k * x(t) = 0

%   [ x1 ]         [   0     1  ]  [ x1 ]
% d [    ] / dt =  [            ]  [    ]
%   [ x2 ]         [ -k/m  -b/m ]  [ x2 ] 

m = 20; % mass
k = 7; % spring constant
b = 2; % damping coefficient
x0 = [1, 0]'; % initial state vector (position and velocity)
R = 0.09; % variance of measurement noise
R = R*100;

dt= 0.25
F = [0 1;-k/m -b/m];
phi = expm(F*dt);
Sigma0=diag([0.1,0.1]);
s=0.0004;
% s=s.*100;
Q = s*diag([dt,dt]);

xnn=x0;
Pnn=Sigma0;
H = [1 0];
I = eye(size(Pnn));
x_KF=[];
P_KF=[];
% figure 
% hold on
for n=1:length(t)
   
    xnn_p= phi*xnn;
    Pnn_p= phi*Pnn*phi'+Q;
    K= Pnn_p*H'*(H*Pnn_p*H'+ R)^(-1);
    xnn=xnn_p+K*(z(n)-xnn_p(1));
    Pnn=(I-K*H)*Pnn_p;
    % save result 
    x_KF=[x_KF xnn];
    Sigma_KF(:,n)=sqrt(diag(Pnn));
    P_KF=[P_KF Pnn];
%     plot(n,diag(Pnn_p)','o')
    
end

%% True position
rhs = @(t,x) F*x; % rhs of function
[~,trueTrajectory] = ode45(rhs,t,[1,0]);

figure
hold on
plot(t,trueTrajectory(:,1));
plot(t,x_KF(1,:));
legend('True Position','Filtered Position')
% plot(t,z,'c');
% legend('True Position','Filtered Podition ','Observation')
title(['Position with Measurement Noise:\sigma^{2}_{0} =  ',num2str(R),'  Process noise:\sigma^{2}_{0} =  ',num2str(s) ])

figure
hold on
plot(t,trueTrajectory(:,2));
plot(t,x_KF(2,:));
legend('True Velocoty','Filtered Velocity')
title(['Velosity with Measurement Noise:\sigma^{2}_{0} =  ',num2str(R),'  Process noise:\sigma^{2}_{0} =  ',num2str(s) ])


figure
subplot(2,1,1)
hold on
plot(t,Sigma_KF(1,:));
title(['Standard Deviation of Position with Measurement Noise:\sigma^{2}_{0} =  ',num2str(R),'  Process noise:\sigma^{2}_{0} =  ',num2str(s) ])
subplot(2,1,2)
hold on
plot(t,Sigma_KF(2,:));
title(['Standard Deviation of Velocity with Measurement Noise:\sigma^{2}_{0} =  ',num2str(R),'  Process noise:\sigma^{2}_{0} =  ',num2str(s) ])