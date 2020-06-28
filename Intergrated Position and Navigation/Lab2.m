% Integrated Positioning and Navigation Lab2
% 3371477
% Tianqi Xiao

% Task 1
%% Rotation angles
W1=deg2rad([0.3 0.2 0.05]);
W2=deg2rad([-30 35 -20]);

%% DCM of equation 2.4
Cst1_1=[ cos(W1(3))*cos(W1(1))-cos(W1(2))*sin(W1(3))*sin(W1(1)), cos(W1(3))*sin(W1(1))+ cos(W1(2))*sin(W1(3))*cos(W1(1)), sin(W1(2))*sin(W1(3));
        -sin(W1(3))*cos(W1(1))-cos(W1(2))*cos(W1(3))*sin(W1(1)), -sin(W1(3))*sin(W1(1))+ cos(W1(2))*cos(W1(3))*cos(W1(1)), sin(W1(2))*cos(W1(3));
         sin(W1(2))*sin(W1(1)),                                  -sin(W1(2))*cos(W1(1)),                                   cos(W1(2))]

Cst1_2=[ cos(W2(3))*cos(W2(1))-cos(W2(2))*sin(W2(3))*sin(W2(1)), cos(W2(3))*sin(W2(1))+ cos(W2(2))*sin(W2(3))*cos(W2(1)), sin(W2(2))*sin(W2(3));
        -sin(W2(3))*cos(W2(1))-cos(W2(2))*cos(W2(3))*sin(W2(1)),-sin(W2(3))*sin(W2(1))+ cos(W2(2))*cos(W2(3))*cos(W2(1)), sin(W2(2))*cos(W2(3));
         sin(W2(2))*sin(W2(1)),                                 -sin(W2(2))*cos(W2(1)),                                   cos(W2(2))]

%% DCM of equation 2.6

Cst2_1=[ cos(W1(2))*cos(W1(3)),                                   cos(W1(2))*sin(W1(3)),                                 -sin(W1(2));
        -sin(W1(3))*cos(W1(1))+sin(W1(2))*cos(W1(3))*sin(W1(1)),  cos(W1(3))*cos(W1(1))+sin(W1(2))*sin(W1(3))*sin(W1(1)), cos(W1(2))*sin(W1(1));
         sin(W1(3))*sin(W1(1))+sin(W1(2))*cos(W1(3))*cos(W1(1)), -cos(W1(3))*sin(W1(1))+sin(W1(2))*sin(W1(3))*cos(W1(1)), cos(W1(2))*cos(W1(1))]

Cst2_2=[ cos(W2(2))*cos(W2(3)),                                    cos(W2(2))*sin(W2(3)),                                  -sin(W2(2));
        -sin(W2(3))*cos(W2(1))+ sin(W2(2))*cos(W2(3))*sin(W2(1)),  cos(W2(3))*cos(W2(1))+ sin(W2(2))*sin(W2(3))*sin(W2(1)), cos(W2(2))*sin(W2(1));
         sin(W2(3))*sin(W2(1))+ sin(W2(2))*cos(W2(3))*cos(W2(1)), -cos(W2(3))*sin(W2(1))+ sin(W2(2))*sin(W2(3))*cos(W2(1)), cos(W2(2))*cos(W2(1))]
%% Euler Symmertic Parameters
% Accroding to (2.11), trace(Cst£©=1+2*cos(phi) f=(diag(Cst)-cos(phi))/(1-cos(phi))

cos_phi=0.5*[trace(Cst1_1)-1 trace(Cst1_2)-1 trace(Cst2_1)-1 trace(Cst2_2)-1];
phi=acos(cos_phi)

f=sqrt(([diag(Cst1_1) diag(Cst1_2) diag(Cst2_1) diag(Cst2_2)]-cos_phi)./(1-cos_phi))

f1=sqrt(([Cst1_1(1,1) Cst1_2(1,1) Cst2_1(1,1)  Cst2_2(1,1)]-cos_phi)./(1-cos_phi));
f2=sqrt(([Cst1_1(2,2) Cst1_2(2,2)  Cst2_1(2,2)  Cst2_2(2,2) ]-cos_phi)./(1-cos_phi));
f3=sqrt(([Cst1_1(3,3) Cst1_2(3,3)  Cst2_1(3,3)  Cst2_2(3,3) ]-cos_phi)./(1-cos_phi));



