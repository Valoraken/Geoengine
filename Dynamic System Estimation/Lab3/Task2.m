clear
% clc
%[x' y']

F_R=[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];

G_R=[1 0;0 0;0 1;0 0];
% [x' x" y' y"]'
F_IR=[0 1 0 0;0 0 0 0;0 0 0 1;0 0 0 0];
G_IR=[0 0;1 0;0 0;0 1];
W_R=[0.64 0;0 0.04];
W_IR=[0.64 0;0 0.04];

Sigma_xy=[1.2 0 0.3 0; 0 0.64  0 0 ; 0.3 0 1.0 0;0 0 0 0.04];
dt=1;

A_R=[-F_R G_R*W_R*G_R'; zeros(size(F_R)) F_R'];%(8 by 8)
A_IR=[-F_IR G_IR*W_IR*G_IR'; zeros(4,4) F_IR'];%(8 by 8)

B_R=expm(A_R);
B_IR=expm(A_IR);

phi_R=B_R((5:8),(5:8))'
phi_IR=B_IR((5:8),(5:8))';

Q_R=phi_R*B_R((1:4),(5:8));
Q_IR=phi_IR*B_IR((1:4),(5:8));
Sigma_xy_t=phi_R*Sigma_xy*phi_R'+Q_R
Sigma_xy_2t=phi_IR*Sigma_xy_t*phi_IR'+ Q_IR 
K=sqrt((Sigma_xy_2t(1,1)-Sigma_xy_2t(3,3))^2+4*Sigma_xy_2t(1,3)^2);
a=sqrt(0.5*(Sigma_xy_2t(1,1)+Sigma_xy_2t(3,3)+K))
b=sqrt(0.5*(Sigma_xy_2t(1,1)+Sigma_xy_2t(3,3)-K))
phi=atan2(2*Sigma_xy_2t(1,3),(Sigma_xy_2t(1,1)-Sigma_xy_2t(3,3)))/2

% phi=expm(F_R*dt)
%% without velosity in RW
F_R0=zeros(2,2);
G_R0=[1 0;0 1];
Sigma_xy0=[1.2 0.3; 0.3 1];
A_R0=[-F_R0 G_R0*W_R*G_R0'; zeros(size(F_R0)) F_R0'];
B_R0=expm(A_R0);
phi_R0=B_R0((3:4),(3:4))'
Q_R0=phi_R0*B_R0((1:2),(3:4));

Sigma_xy_dt=phi_R0*Sigma_xy0*phi_R0'+Q_R0
Sigma_xy_dt2=[Sigma_xy_dt(1,1) 0 Sigma_xy_dt(1,2) 0; 0 W_IR(1,1) 0 0;Sigma_xy_dt(2,1) 0 Sigma_xy_dt(2,2) 0;0 0 0 W_IR(2,2)];
Sigma_xy_2=phi_IR*Sigma_xy_dt2*phi_IR'+ Q_IR

K0=sqrt((Sigma_xy_2(1,1)-Sigma_xy_2(3,3))^2+4*Sigma_xy_2(1,3)^2);
a0=sqrt(0.5*(Sigma_xy_2(1,1)+Sigma_xy_2(3,3)+K))
b0=sqrt(0.5*(Sigma_xy_2(1,1)+Sigma_xy_2(3,3)-K))
phi0=atan2(2*Sigma_xy_2(1,3),(Sigma_xy_2(1,1)-Sigma_xy_2(3,3)))/2
