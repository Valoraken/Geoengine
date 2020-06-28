clear
close all
% random walk
F=0;
dt=1;
phi=expm(F*dt);
Q=4*dt;
t=1:200;
Sigma0=100;

R=1;
x1=0;
x2=0;
I=1;
P1=Sigma0;
P2=Sigma0;
P1_KF=[];
P2_KF=[];
for n=1:length(t)
 
%     x1_p= phi*x1;
    P1_p= phi*P1*phi'+Q;
    H1 = 0.5;
    K1= P1_p*H1'*(H1*P1_p*H1'+ R)^(-1);
%     x1=x1_p+K1*(z1-H1*x1_p);
    P1=(I-K1*H1)*P1_p;
    P1_KF=[P1_KF P1];
    Sigma_1(n)=sqrt(diag(P1));

    
%     x2_p= phi*x2;
    P2_p= phi*P2*phi'+Q;
    H2 =  cos(1+n/90);
    K2= P2_p*H2'*(H2*P2_p*H2'+ R)^(-1);
%     x2=x2_p+K2*(z2-sqrt(x2_p^2+h^2));
    P2=(I-K2*H2)*P2_p;
    P2_KF=[P2_KF P2];
    Sigma_2(n)=sqrt(diag(P2));
    K_2(n)=K2;
    KH_2(n)=K2*H2;
end
figure
hold on
plot(t,Sigma_1);
plot(t,Sigma_2);
legend('Instrument 1 ','Instrument 2 ')
figure
hold on
plot(t,meshgrid(H1,t)');
plot(t,cos(1+t./90));
legend('H1','H2')
figure
hold on
plot(t,meshgrid(K1,t)');
plot(t,K_2);
legend('K1','K2')