% DSE Lab2
% Task 1
close all
clear all
clc


dt=0.1;
t=0:dt:10;

% generate Gaussian white noise
W=normrnd(0,1,1,max(t)/dt);

%% case: x"=0+W(t)
Xi=[0 0]';
F_0=[0,1;0,0];
Phi_0=expm(F_0*dt);
G=[0;1];
sigma=1;
X=zeros(2,(max(t)/dt)+1);
X(:,1)=Xi;
Sigma_x=zeros(1,(max(t)/dt)+1);
Sigma_xxi=zeros(2,2);
Sigma_x(1)=sqrt(Sigma_xxi(1,1));
Q=[dt^3/3 dt^2/2;dt^2/2 dt];
for i=1:max(t)/dt
    Xi=Phi_0*Xi+[0;W(i)];
    X(:,i+1)=Xi;
    Sigma_xxi=Phi_0*Sigma_xxi*Phi_0'+Q;
%     Sigma_x(i+1)=sqrt(Sigma_xxi(1,1));
    Sigma_x(i+1)=(Sigma_xxi(1,1));

end

%% case: xc"=0+c(t) with c'=-beta*c(t)+W(t)
% Define beta(0,1) and transition function phi
beta = 0.7;
% beta1=input('Input a number between(0,1) as beta = ')
% transition function
F=[0 1;0 0]
b=[0;1];
Fc=-beta;
F_new=[F b;zeros(size(b')) Fc];
Phi_1=expm(F_new*dt);
G=[0 0;0 0];
Gc=1;
G_new=[G zeros(size (b));zeros(size(b')) Gc];
Xci=[0;0;0];
Xc=zeros(3,(max(t)/dt)+1);
Xc(:,1)=Xci;
Sigma_xc=zeros(1,(max(t)/dt)+1);
Sigma_xcxci=zeros(3,3);
Sigma_xc(1)=sqrt(Sigma_xcxci(1,1));
W_new=[zeros(size(W));zeros(size(W));W];
GWGT=G_new*G_new';
A=[-F_new GWGT;zeros(size(GWGT)) F_new']*dt;
B=expm(A);
Qc=B(4:6,4:6)'*B(1:3,4:6);
for i=1:max(t)/dt
    Xci=Phi_1*Xci+Phi_1*G_new*W_new(:,i);
    Xc(:,i+1)=Xci;
    Sigma_xcxci=Phi_1*Sigma_xcxci*Phi_1'+Qc;
%     Sigma_xc(i+1)=sqrt(Sigma_xcxci(1,1));
    Sigma_xc(i+1)=(Sigma_xcxci(1,1));
end

%% visualize
figure (1)
subplot(2,1,1)
hold on
% scatter(t(1:max(t)/dt),W);
plot(t,X(1,:));
plot(t,X(2,:));
% legend('w ','x ','v')
legend('x(t)','v')
subplot(2,1,2)
hold on
% scatter(t(1:max(t)/dt),W);
% plot(t,Xc(3,:));
plot(t,Xc(1,:));
plot(t,Xc(2,:));
% legend('w ',' c','x ','v')
legend('xc(t) ','v');
figure(2)
hold on
plot(t,Sigma_x);
plot(t,Sigma_xc);
legend('varience - x(t)','varience - xc(t)')
