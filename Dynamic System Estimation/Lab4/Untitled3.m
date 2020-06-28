% DSE Lab4
% Xiao Tianqi
% 3371477
% 17.6.2019
clear
close all
%% Observations
Z=load('exampleKF2.txt')';
% x Vx
Z=Z(1:2,:);
H=[1 0 0 0;
   0 1 0 0];
% covariance
R=eye(2,2)*50*50;
%% Prediction
n = 30;
GM=3.986005*1e14;
x0=[-12823317 -11933101 2000 -1000]';
xnn=x0;
G=n*[0 0 1 1]';
Sigma0=diag([10,10,100,100]);
Pnn=Sigma0;
x_pre=[];
x_fwd=[];
P_pre=[];
P_fwd=[];
dt = 1;
I = eye(4,4);
%% Forward Karman filter
for i = 1:length(Z)
a=-GM/(norm(xnn(1:2))^3);
F=[ 0 0 1 0;
    0 0 0 1;
    a 0 0 0;
    0 a 0 0];
% get Q and Phi
A=dt*[-F G*G';
    zeros(size(F)) F'];
B=expm(A);
Phi = B(1+length(A)/2:length(A),1+length(A)/2:length(A))';
Q= Phi*B(1:length(A)/2,1+length(A)/2:length(A));
% apply filter
    xnn_p= Phi*xnn;
    Pnn_p= Phi*Pnn*Phi'+Q;
    K= Pnn_p*H'/(H*Pnn_p*H'+ R);
    xnn=xnn_p+K*(Z(:,i)-H*xnn_p);
    Pnn=(I-K*H)*Pnn_p;
    % save result 
    x_pre=[x_pre xnn_p];
    x_fwd=[x_fwd xnn];
    Sigma_fwd(:,i)=sqrt(diag(Pnn));
    P_fwd=[P_fwd Pnn];
 
end

%% Backward Karman Filter
xnn =x_fwd(:,length(x_fwd));
Pnn =eye(4,4); 
x_bwd=[];
P_bwd=[];
for i = 1:length(Z)
a=-GM/(norm(xnn(1:2))^3);
F=[ 0 0 1 0;
    0 0 0 1;
    a 0 0 0;
    0 a 0 0];
% get Q and Phi
A=[-F G*G';zeros(size(F)) F'];
B=expm(A);
Phi = B(1+length(A)/2:length(A),1+length(A)/2:length(A))';
Q= Phi*B(1:length(A)/2,1+length(A)/2:length(A));
% apply filter
    xnn_p= Phi*xnn;
    Pnn_p= Phi*Pnn*Phi'+Q;
    K= Pnn_p*H'/(H*Pnn_p*H'+ R);
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
    
    Pf= P_fwd(:,(i-1)*4+1:i*4)
    Pb= P_bwd(:,(i-1)*4+1:i*4)
    Pn= (Pf^-1+Pb^-1)^-1
    Xn= Pn*(Pf^-1*x_fwd(:,i)+Pb^-1*x_bwd(:,i))
    x=[x Xn];
    xf=x_fwd(:,i)
    xb=x_bwd(:,i)
    Sigma=[Sigma sqrt(diag(Pn))];
    
end

%% plot
% figure
% hold on
% plot(1:Z(1,:));
% plot(x_fwd(1,:),x_fwd(2,:),x_fwd(3,:),'o');
% plot(x_bwd(1,:),x_bwd(2,:),x_bwd(3,:),'o');
% plot(x(1,:),x(2,:),x(3,:),'o')
% legend('Z','Forward ','Backward','smooth ');
% axis normal
% grid on
% xlabel('X');ylabel('Y');zlabel('Z')

figure
hold on
 plot(1:length(Z),Z(1,:),'o');
plot(1:length(Z),x_fwd(1,:),'o');
plot(1:length(Z),x_bwd(1,:),'*');
plot(1:length(Z),x(1,:));
legend('Z','Forward ','Backward','smooth ');
axis normal
grid on
xlabel('X');ylabel('Y');zlabel('Z')

figure
hold on
plot(1:length(Z),x_fwd(2,:),'*');
plot(1:length(Z),x_bwd(2,:),'o');
plot(1:length(Z),x(2,:));
legend('Forward ','Backward','smooth ');
axis normal
grid on
xlabel('X');ylabel('Y');zlabel('Z')

% figure
% hold on
% plot(1:length(Z),x_fwd(3,:));
% plot(1:length(Z),x_bwd(3,:));
% plot(1:length(Z),x(3,:));