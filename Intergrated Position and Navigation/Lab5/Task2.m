dv=load('task2_accel.txt');
N=1000;
f=100;
dt=2*1/f;

for tk=2:N
    a1(tk)=(3*dv(tk)-dv(tk-1))/dt;
    a2(tk-1)=(dv(tk-1)+dv(tk))/dt;
end
for tk=3:N
    a3(tk-2)=(3*dv(tk-1)-dv(tk))/dt;
end
t=0:1/f:N/f;
plot(t(2:end),a1)
hold on;
plot(t(2:end-1),a2)
plot(t(2:end-2),a3)
legend('a^p(t_k)','a^p(t_k_-_1)','a^p(t_k_-_2)')