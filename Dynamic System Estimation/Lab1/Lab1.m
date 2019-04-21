Y= [4.96 4.93 5.25 -5.26 -4.96 -4.99;4.96 4.94 5.56 -5.55 -4.89 -4.92;4.97 4.91 5.84 -5.81 -4.91 -4.92
    ;4.86 4.83 6.09 -6.09 -4.86 -4.89;4.85 4.81 6.33 -6.34 -4.80 -4.82]'
A=[1 0 0;-1 1 0; 0 -1 1; 0 1 -1;1 -1 0;-1 0 0]
P=eye(6,6)
% separate
for i=1:5
    Hs(:,i)=A\Y(:,i)
end
%sequencial
x_hat=inv(A'*P*A)*A'*P*Y(:,1)
e_1=Y(:,1)-A*x_hat
s_02=e_1'*P*e_1/(6-3)
Sxx=s_02^2*inv(A'*P*A)
sigma(1,1)=Sxx(1,1)
sigma(2,1)=Sxx(2,2)
sigma(3,1)=Sxx(3,3)
e(1,:)=e_1;
for i=2:5
   x_old=x_hat;
   
   s_old=s_02;
   Sxx_old=Sxx;
   yk=Y(:,i);
   x_hat=x_old+inv(s_old*inv(Sxx_old)+A'*P*A)*A'*P*(yk-A*x_old)
   dx=x_hat-x_old;
   s_02=s_old*(6-3+dx'*inv(Sxx_old)*dx)+(yk-A*x_old)'*P*(yk-A*x_old);
   Sxx=s_02*inv(s_old*inv(Sxx_old)+A'*P*A);
   
end
