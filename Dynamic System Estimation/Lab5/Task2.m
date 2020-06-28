clear
close all


% observations
data=load('EFK_task2.txt');
t=data(:,1)';
DP1=data(:,2)';         % distance to P1
DP2=data(:,3)';         % distance to P2
dt = t(2)-t(1);
R= [0.01^2 0 ;0 0.01^2 ]
Z=[DP1;DP2];
% control points
P1=[-3.5;-3.5];
P2=[3.5;-3.5];

% true positions
%%% x y of First Given True position is reversed compared to the given initial value x0
r= 2+sin(20*pi.*t);
X = r.*sin(2*pi.*t);
Y = r.*cos(2*pi.*t);
Vx = 20*pi*cos(20*pi.*t).*sin(2*pi.*t) + 2*pi*cos(2*pi.*t).*(sin(20*pi.*t) + 2);
Vy = 20*pi*cos(2*pi.*t).*cos(20*pi.*t) - 2*pi*sin(2*pi.*t).*(sin(20*pi.*t) + 2);

% intergrated random walk: 
% x" = 0 + w(t)  y" = 0 + w(t)
% x = [x y x' y']'
x0=[2 0 0 0]';
F=[0 0 1 0;
   0 0 0 1;
   0 0 0 0;
   0 0 0 0];
 Phi=expm(F*dt); 
% Phi =[1 dt 0 0;
%       0 1 0 dt;
%       0 0 1  0;
%       0 0 0 1];
Sigma0=diag([10,10,10,10]);
Q = 4*[dt^3/3     0      dt^2/2      0  ;
          0     dt^3/3     0      dt^2/2;
       dt^2/2     0        dt        0  ;
          0     dt^2/2     0         dt];
% DP1 = sqrt((x-P1x)^2+(y-P1y)^2)
% DP2 = sqrt((x-P2x)^2+(y-P2y)^2)
xnn=x0;
Pnn=Sigma0;
x_pre=[];
x_fwd=[];
P_fwd=[];
I = eye(size(Sigma0));
for i=1:length(t)
  xnn_p= Phi*xnn;
  Pnn_p= Phi*Pnn*Phi'+Q;
   
% Jacobian
d=[sqrt((xnn_p(1)-P1(1))^2+(xnn_p(2)-P1(2))^2);
   sqrt((xnn_p(1)-P2(1))^2+(xnn_p(2)-P2(2))^2)];
J = [(xnn_p(1)-P1(1))/d(1) (xnn_p(2)-P2(2))/d(1);
     (xnn_p(1)-P2(1))/d(2) (xnn_p(2)-P2(2))/d(2)];
H = [J zeros(size(J))];
    K= Pnn_p*H'*(H*Pnn_p*H'+ R)^(-1);
    xnn=xnn_p+K*(Z(:,i)-d);
    Pnn=(I-K*H)*Pnn_p;  
    % save result 
    x_pre=[x_pre xnn_p];
    x_fwd=[x_fwd xnn];
    Sigma_fwd(:,i)=sqrt(diag(Pnn));
    P_fwd=[P_fwd Pnn];
end
%% Plot
figure(1)
hold on
plot(X,Y);
plot(x_fwd(1,:),x_fwd(2,:)); 
xlabel('x [m]')
ylabel('y [m]')
legend('True Position (X,Y)','Filtered Position (x,y) ')
title('Position of the  vehicle')
axis equal
grid on



figure(2)
subplot(2,1,1)
hold on
plot(t,X);
plot(t,x_fwd(1,:)); 
xlabel('t [s] ')
ylabel('x [m]')
legend('True Coorsinates X','Filtered Coordinates x ')
title('X coordinates of the vehicle')

subplot(2,1,2)
hold on
plot(t,Y);
plot(t,x_fwd(2,:)); 
xlabel('t [s] '); 
ylabel('y [m]')
legend('True Coorsinates Y','Filtered Coordinates y')
title('Y coordinates of the vehicle')

figure(3)
subplot(2,1,1)
hold on
plot(t,Vx);
plot(t,x_fwd(3,:)); 
xlabel('t [s] ')
ylabel('Vx [m/s]')
legend('True Velocity Vx','Filtered Velocity vx ')
title('Velocity of the vehicle in x direction')

subplot(2,1,2)
hold on
plot(t,Vy);
plot(t,x_fwd(4,:)); 
xlabel('t [s] '); 
ylabel('Vy [m/s]')
legend('True Velocity Vy','Filtered Velocity vy ')
title('Velocity of the vehicle in y direction')


% standard deviation
figure(4)
hold on 
title('Standard deviation of unknowns ')
subplot(1,2,1)
hold on
plot(t, Sigma_fwd(1:2,:)');
legend('\sigma x','\sigma y');
subplot(1,2,2)
hold on
plot(t, Sigma_fwd(3:4,:)');
legend('\sigma vx','\sigma vy');
