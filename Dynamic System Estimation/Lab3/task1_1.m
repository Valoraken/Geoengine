clc
clear
B=0.7;
F=[0 1 0;0 0 1;0 0 -B];
G=[0;0;1];
Z=zeros(3,3);
dt=0.1;
W=[1];
Q1=zeros(3,3);
Q2=zeros(1,100);
for i=1:100
A=[-F G*W*G';Z F'].*dt;
B=expm(A);
phi=B(4:6,4:6)';
Q=phi*B(1:3,4:6);
Q1=phi*Q1*phi'+Q;
Q2(i)=Q1(1,1);
end
%-------------------
Q3=zeros(2,2);
Q4=zeros(1,100);
F2=[0 1;0 0];
G2=[0;1];
Z2=zeros(2,2);
for j=1:100
A2=[-F2 G2*W*G2';Z2 F2'].*dt;
B2=expm(A2);
phi=B2(3:4,3:4)';
q=phi*B2(1:2,3:4);
Q3=phi*Q3*phi'+q;
Q4(j)=Q3(1,1);
end
plot(Q2,'r-');hold on
plot(Q4,'g-');
