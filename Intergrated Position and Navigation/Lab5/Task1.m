S1=load('task1_sensor1.txt'); 
S2=load('task1_sensor2.txt');
N1=length(S1);
N2=length(S2);
f1=100;
f2=150;
[T1,sigma1] = allan(S1,f1,1000);
[T2,sigma2] = allan(S2,f2,1500);
% load('task_res.mat')
figure()
loglog(T1,sigma1);
hold on 
loglog(T2,sigma2)
legend('sensor 1','sensor 2');
xlabel('\tau [s]');
ylabel('\sigma_y');