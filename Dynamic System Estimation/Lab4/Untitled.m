clear
% initial values
Ob=load('exampleKF2.txt');
X_O=Ob(:,1);
Y_O=Ob(:,2);
Z_O=Ob(:,3);
GM=3.986005e14;
x0=[-12823317;-11933101;20070042;2000;-1000;1000];
W=1;
G=30*[0;0;0;1;1;1];
x1=[x0 zeros(6,501)];
P1=[x0 zeros(6,501)];
Rn=2500*eye(3,3);
H=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0];
P=[10 0 0 0 0 0;0 10 0 0 0 0;0 0 10 0 0 0;0 0 0 100 0 0;0 0 0 0 100 0; 0 0 0 0 0 100];
I=eye(6,6);
% a=zeros(501);
for i=1:501
tx=x0(1);ty=x0(2);tz=x0(3);
a=-3.986005e14/(tx^2+ty^2+tz^2)^(3/2);
F=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1;a 0 0 0 0 0;0 a 0 0 0 0;0 0 a 0 0 0];
phi=expm(F);
A1=[-F G*W*G';zeros(6,6) F'];
B1=expm(A1);
B22=B1(7:12,7:12);
B12=B1(1:6,7:12);
% phi=B22';
Q1=phi*B12;%
x0=phi*x0;
P=phi*P*phi'+Q1;
% Z(:,i)=H*x1(:,i+1)+[50;50;50];
Z=Ob';
% K=P*H'*inv(H*P*H'+Rn);
K=P*H'*inv(H*P*H'+Rn);
% x0=x0+K*(Z(:,i)-H*x0);
% P=(I-K*H)*P;
% x2(:,i)=x1(:,i+1)+K*(Z-H*x1(:,i+1));
% P2=()
P=(I-K*H)*P;
P_for{i}=P;
x0=x0+K*(Ob(i,:)'-H*x0);
x_forward=x0;
x_for{i}=x_forward;
X_forward(i)=x_forward(1);
Y_forward(i)=x_forward(2);
Z_forward(i)=x_forward(3);
end
%% backward kalman filter
dt= 1;
for i=1:501
tx=x0(1);ty=x0(2);tz=x0(3);
a=-3.986005e14/(tx^2+ty^2+tz^2)^(3/2);
F=[0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;a,0,0,0,0,0;0,a,0,0,0,0;0,0,a,0,0,0];
G=30*[0;0;0;1;1;1];
phi=expm(F*dt);

A1=[-F,G*W*G';zeros(6,6),F']*dt;
B1=expm(A1);
%phi=B(7:12,7:12)';
Q1=phi*B1(1:6,7:12);

% prediction
x0=phi*x0;
P=phi*P*phi'+Q1;
% kalman gain
K=P*H'*(H*P*H'+Rn)^(-1);
% update
I=[1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1];
P=(I-K*H)*P;
P_back{501-i+1}=P;
x0=x0+K*(Ob(501-i+1,:)'-H*x0);
x_backward=x0;
x_back{501-i+1}=x_backward;
X_backward(501-i+1)=x_backward(1);
Y_backward(501-i+1)=x_backward(2);
Z_backward(501-i+1)=x_backward(3);
end

%% smoothed kalman filter
for i=1:501
    AA=P_back{i}*(P_for{i}+P_back{i})^(-1);
    x_smooth=AA*x_for{i}+(I-AA)*x_back{i};
    P_smooth=AA*P_for{i}*AA'+(I-AA)*P_back{i}*(I-AA)';
    X_smooth(i)=x_smooth(1);
    Y_smooth(i)=x_smooth(2);
    Z_smooth(i)=x_smooth(3);
end

%% plot result
figure(5)
plot3(X_O,Y_O,Z_O,'k'); hold on;
plot3(X_forward,Y_forward,Z_forward,'r'); hold on;
plot3(X_backward,Y_backward,Z_backward,'g'); hold on;
plot3(X_smooth,Y_smooth,Z_smooth,'b'); 
legend('observation','forward','backward','smooth');