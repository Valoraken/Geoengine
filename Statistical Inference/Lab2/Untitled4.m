hm=[6.9803+77/1000,10.0213-77/1000,14.3030,4.2871,7.5001+77/1000]'
HA=100.956;
HB=100.459;
% case 1: HA is given, the number of unknowns are 4
% reduced obsevaion y 
y1=hm;
y1(5)=y1(5)-HA 
% design matrix A
A1=[-1 0 0 1;-1 0 1 0;-1 1 0 0;0 1 -1 0;-1 0 0 0]
rank(A1)
Ay1=[A1 y1]
rank(Ay1)
% N1=A1'*A1
% x1hat=inv(N1)*A1'*y1
% Least-suqre estimate of the vector Hhat
x1hat=A1\y1
% corresponding adjusted reduced observation yhat
y1hat=A1*x1hat
% adjusted observations hhat
h1hat=y1hat
h1hat(5)=h1hat(5)+HA
% edstimate ehat
e1hat=y1-y1hat
% quardic sum ehat'ehat
sum1=e1hat'*e1hat
% redundancy r
r1=1
sumr1=sum1/r1
%orth check
or1=A1'*e1hat
%main check
mainc1=y1hat-A1*x1hat
%case 2

y2=hm;
y2(1)=y2(1)-HB 
% design matrix A
A2=[0 -1 0 0;0 -1 0 1 ;0 -1 1 0 ;0 0 1 -1 ;1 -1 0 0]
rank(A2)
Ay2=[A2 y2]
rank(Ay2)
% N2=A2'*A2
% x2hat=inv(N2)*A2'*y2
% Least-suqre estimate of the vector Hhat
x2hat=A2\y2
% corresponding adjusted reduced observation yhat
y2hat=A2*x2hat
% adjusted observations hhat
h2hat=y2hat
h2hat(1)=h2hat(1)+HB
% edstimate ehat
e2hat=y2-y2hat
% quardic sum ehat'ehat
sum2=e2hat'*e2hat
% redundancy r
r2=1
sumr2=sum2/r2
%orth check
or2=A2'*e2hat
%main check
mainc2=y2hat-A2*x2hat
%case 3
A3=[-1 0 0 ;-1 0 1 ;-1 1 0 ;0 1 -1 ;-1 0 0 ]
rank(A3)
y3=hm;
y3(5)=y3(5)-HA;y3(1)=y3(1)-HB
Ay3=[A3 y3]
rank(Ay3)
% N3=A3'*A3
% x3hat=inv(N3)*A3'*y3
x3hat=A3\y3
y3hat=A3*x3hat
% adjusted observations hhat
h3hat=y3hat
h3hat(1)=h3hat(1)+HB
h3hat(5)=h3hat(5)+HA
% edstimate ehat
e3hat=y3-y3hat
% quardic sum ehat'ehat
sum3=e3hat'*e3hat
% redundancy r
r3=2
sumr3=sum3/r3
%orth check
or3=A3'*e3hat
%main check
mainc3=y3hat-A3*x3hat