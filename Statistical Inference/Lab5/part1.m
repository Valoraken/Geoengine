k=77;

d=[100.04+k/1000 100.1 99.98 200.0 200.02 299.96+k/1000]'; 
A=[1 0 0;0 1 0;0 0 1;1 1 0;0 1 1;1 1 1];
B=[1 1 0 -1 0 0;0 1 1 0 -1 0;1 1 1 0 0 -1]';
% orthogonlal check: B'*A
%% compute A model
x_hatA=inv(A'*A)*A'*d
d_hatA=A*x_hatA
e_hatA=d-d_hatA

%% compute B model

e_hatB=B*inv(B'*B)*B'*d
d_hatB=d-e_hatB