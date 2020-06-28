TE=365.356;% d
TM=686.971;% d
aE=149.5980e6;% km
aM=227.9392e6;% km
(TM^2/aM^3)-(TE^2/aE^3)
y=[TM;aE;TE;aM];
e0=[0;0;0;0];%eTM;eaE;eTE;eaM
delta =1;
while delta>10e-10
fe0=(TM-e0(1))^2*(aE-e0(2))^3-(TE-e0(3))^2*(aM-e0(4))^3
B=[];
w=fe0-B'*e0;
e_hat=B*inv(B'*B)*w;
delta=e_hat-e0;
e0=e_hat;
end