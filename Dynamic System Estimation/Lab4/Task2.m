% DSE Lab4
% Xiao Tianqi
% 3371477
% 17.6.2019

close all
clear
clc
%% Observations
Z=load('exampleKF2.txt')';

H=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0];
% covariance
R=eye(3,3)*50*50;
%% Prediction
n = 30;
GM=3.986005*1e14;
x0=[-12823317 -11933101 20070042 2000 -1000 1000]';
xnn=x0;
G=n*[0 0 0 1 1 1]';
Sigma0=diag([10,10,10,100,100,100]);
Pnn=Sigma0;
x_pre=[];
x_fwd=[];
P_pre=[];
P_fwd=[];
I = eye(6,6);
%% Forward Karman filter
dt = 1;
for i = 1:length(Z)
a=-GM/(sqrt(xnn(1)^2+xnn(2)^2+xnn(3)^2)^3);
F=[ 0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1
    a 0 0 0 0 0
    0 a 0 0 0 0
    0 0 a 0 0 0 ];
% get Q and Phi
A=[-F G*G';
    zeros(size(F)) F']*dt;
B=expm(A);
Phi = B(1+length(A)/2:length(A),1+length(A)/2:length(A))';
Q= Phi*B(1:length(A)/2,1+length(A)/2:length(A));
% apply filter
    xnn_p= Phi*xnn;
    Pnn_p= Phi*Pnn*Phi'+Q;
    K= Pnn_p*H'*(H*Pnn_p*H'+ R)^(-1);
    xnn=xnn_p+K*(Z(:,i)-H*xnn_p);
    Pnn=(I-K*H)*Pnn_p;
    % save result 
    x_pre=[x_pre xnn_p];
    x_fwd=[x_fwd xnn];
    Sigma_fwd(:,i)=sqrt(diag(Pnn));
    P_fwd=[P_fwd Pnn];
 
end

%% Backward Karman Filter
dt = 1;
% xnn =x_fwd(:,length(x_fwd));
% Pnn =Sigma0; 
x_bwd=[];
P_bwd=[];
for i = 1:length(Z)
a=-GM/(sqrt(xnn(1)^2+xnn(2)^2+xnn(3)^2)^3);
F=[ 0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1
    a 0 0 0 0 0
    0 a 0 0 0 0
    0 0 a 0 0 0 ];
% get Q and Phi
A=[-F G*G';zeros(size(F)) F']*dt;
B=expm(A);
Phi = B(1+length(A)/2:length(A),1+length(A)/2:length(A))';% B22
Q= Phi*B(1:length(A)/2,1+length(A)/2:length(A)); % B12
% apply filter
    xnn_p= Phi*xnn;
    Pnn_p= Phi*Pnn*Phi'+Q;
    K= Pnn_p*H'*(H*Pnn_p*H'+ R)^-1;
    xnn=xnn_p+K*(Z(:,length(Z)-i+1)-H*xnn_p);
    Pnn=(I-K*H)*Pnn_p;
    x_bwd=[xnn x_bwd];
    P_bwd=[Pnn P_bwd];
    Sigma_bwd(:,length(Z)-i+1)=sqrt(diag(Pnn));

end

%% smooth filter
x=[];
Sigma=[];
for i = 1:length(Z)
    Pf= P_fwd(:,(i-1)*6+1:i*6);
    Pb= P_bwd(:,(i-1)*6+1:i*6);
    Pn= (Pf^-1+Pb^-1)^-1;
    Xn= Pn*(Pf^-1*x_fwd(:,i)+Pb^-1*x_bwd(:,i));
    x=[x Xn];
%     xf=x_fwd(:,i)
%     xb=x_bwd(:,i)
    Sigma=[Sigma sqrt(diag(Pn))];
    
end

%% Plot Coordinates
% 3D View
figure
hold on
plot3(Z(1,:)',Z(2,:)',Z(3,:)');
plot3(x_fwd(1,:)',x_fwd(2,:)',x_fwd(3,:)');
plot3(x_bwd(1,:)',x_bwd(2,:)',x_bwd(3,:)');
plot3(x(1,:)',x(2,:)',x(3,:)')
legend('Z','Forward ','Backward','smooth ');
axis normal
grid on
% view(37.5, 30)
xlabel('X');ylabel('Y');zlabel('Z')
title('Locations')
% Each direction
figure
subplot(3,1,1)
hold on
plot(1:length(Z),x_fwd(1,:));
plot(1:length(Z),x_bwd(1,:));
plot(1:length(Z),x(1,:));
xlabel('t');ylabel('X')
legend('Forward ','Backward','smooth ');
title('X Coordinates')

subplot(3,1,2)
hold on
plot(1:length(Z),x_fwd(2,:));
plot(1:length(Z),x_bwd(2,:));
plot(1:length(Z),x(2,:));
xlabel('t');ylabel('Y')
legend('Forward ','Backward','smooth ');
title('Y Coordinates')

subplot(3,1,3)
hold on
plot(1:length(Z),x_fwd(3,:));
plot(1:length(Z),x_bwd(3,:));
plot(1:length(Z),x(3,:));
xlabel('t');ylabel('Z')
legend('Forward ','Backward','smooth ');
title('Z Coordinates')

% Plot Velocity
figure
hold on

subplot(3,1,1)
hold on
plot(1:length(Z),x_fwd(4,:));
plot(1:length(Z),x_bwd(4,:));
plot(1:length(Z),x(4,:));
title('Vx ')
xlabel('t');ylabel('Vx')
legend('Forward ','Backward','smooth ');

subplot(3,1,2)
hold on
plot(1:length(Z),x_fwd(5,:));
plot(1:length(Z),x_bwd(5,:));
plot(1:length(Z),x(5,:));
legend('Forward ','Backward','smooth ');
xlabel('t');ylabel('Vy')
title('Vy ')

subplot(3,1,3)
hold on
plot(1:length(Z),x_fwd(6,:));
plot(1:length(Z),x_bwd(6,:));
plot(1:length(Z),x(6,:));
legend('Forward ','Backward','smooth ');
xlabel('t');ylabel('Vz')
title('Vz ')

% Plot variance
figure
subplot(3,1,1)
hold on
plot(1:length(Z),Sigma_fwd(1,:));
plot(1:length(Z),Sigma_bwd(1,:));
plot(1:length(Z),Sigma(1,:));
legend('\sigma x(Forward) ','\sigma x(Backward)','\sigma x(smooth) ');

subplot(3,1,2)
hold on
plot(1:length(Z),Sigma_fwd(2,:));
plot(1:length(Z),Sigma_bwd(2,:));
plot(1:length(Z),Sigma(2,:));
legend('\sigma y(Forward) ','\sigma y(Backward)','\sigma y(smooth) ');

subplot(3,1,3)
hold on
plot(1:length(Z),Sigma_fwd(3,:));
plot(1:length(Z),Sigma_bwd(3,:));
plot(1:length(Z),Sigma(3,:));
legend('\sigma z(Forward) ','\sigma z(Backward)','\sigma z(smooth) ');
