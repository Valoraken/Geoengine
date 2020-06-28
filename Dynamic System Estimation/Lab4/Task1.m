% DSE Lab4
% Xiao Tianqi
% 3371477
% 17.6.2019
clear
close all
n=5

%% Intergrated radom walk Model

F =[ 0 1 0 0 0 0 
     0 0 0 0 0 0
     0 0 0 1 0 0
     0 0 0 0 0 0
     0 0 0 0 0 1
     0 0 0 0 0 0 ]; 
Phi = expm(F);
X0=[5 0 10 0 15 0]';
Sigma0=eye(6,6);
Q= 0.01*[1/3  1/2   0    0   0   0;
         1/2   1    0    0   0   0;
          0    0   1/3  1/2  0   0;  
          0    0   1/2   1   0   0;
          0    0    0    0  1/3 1/2;
          0    0    0    0  1/2  1];

%% Observation equation
Z =[ 4.96 4.93 5.25 -5.26 -4.96 -4.99;
     4.96 4.94 5.56 -5.55 -4.89 -4.92;
     4.97 4.91 5.84 -5.81 -4.91 -4.92;
     4.86 4.83 6.09 -6.09 -4.86 -4.89;
     4.85 4.81 6.33 -6.34 -4.80 -4.82 ]'
  
H =[ 1  0  0  0  0  0
    -1  0  1  0  0  0
     0  0 -1  0  1  0
     0  0  1  0 -1  0
     1  0 -1  0  0  0
    -1  0  0  0  0  0 ];

% H=[1 0 0;-1 1 0; 0 -1 1; 0 1 -1;1 -1 0;-1 0 0]';
R =0.02*0.02.*eye(6,6);
Xnn=X0;
Pnn=Sigma0;
I=eye(6,6);
X_p=[];
X_cor=[];
P_p=[];
P_cor=[];
%% Karman Filter
for i = 1:n 
    Xnn_p= Phi*Xnn;
    Pnn_p= Phi*Pnn*Phi'+Q;
    Kn= Pnn_p*H'*inv(H*Pnn_p*H'+ R);
    Xnn=Xnn_p+Kn*(Z(:,i)-H*Xnn_p);
    Pnn=(I-Kn*H)*Pnn_p;
    X_p=[X_p Xnn_p];
    X_cor=[X_cor Xnn];
    Sigma(:,i)=sqrt(diag(Pnn));
    P_cor=[P_cor Pnn];
    Xnn =Xnn;
    Pnn =Pnn;
end
%% Location of each epoch
H_pre=X_p(1:2:6,:)';
H_cor=X_cor(1:2:6,:)';
figure
hold on
plot(H_pre);
% legend('Predicted H1','Predicted H2','Predicted H3');
plot(H_cor);
legend('Predicted H1','Predicted H2','Predicted H3','Filtered H1','Filtered H2','Filtered H3');

V_pre=X_p(2:2:6,:)';
V_cor=X_cor(2:2:6,:)';

figure
subplot(1,2,1)
hold on
plot(V_pre);
legend('Predicted V1','Predicted V2','Predicted V3');
subplot(1,2,2)
hold on
plot(V_cor,'-+');
legend('Filtered V1','Filtered V2','Filtered V3');

%% varience
figure
subplot(1,2,1)
hold on
plot(Sigma(1:2:6,:)');
legend('\sigma H1 ','\sigma H2 ','\sigma H3 ')
subplot(1,2,2)
hold on
plot(Sigma(2:2:6,:)');
legend('\sigma v1 ','\sigma v2 ','\sigma v3 ')