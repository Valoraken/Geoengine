% DSE Lab2
% Task 1
close all
clear all
clc

t=1:100;
dt=1;
x0=0;

%% Define beta(0,1) and transition function phi
beta0 = 0.001;
beta1=input('Input a number between(0,1) as beta = ')
% transition function
phi0=exp(-beta0*dt);
phi1=exp(-beta1*dt);

%% realization
for j=1:30
    % generate random number
    W=normrnd(0,1,1,100);
    % initialize temp variable
    xi_1=x0;
    xi_2=x0;
    %% time steps
    for i=0:99
        
        x_1(j,i+1)=phi0*xi_1+W(i+1);
        x_2(j,i+1)=phi1*xi_2+W(i+1);
        
        xi_1=x_1(j,i+1);
        xi_2=x_2(j,i+1);
        
    end
    beta2=['beta =',num2str(beta1)];
%     %    plot result
%     figure(ceil(j/10));
%     subplot(2,5,(mod(j-1,10)+1));
%     hold on;
%     plot(t,x_2(j,:));
%     plot(t,x_1(j,:));
%     plot(t,W,'.');
%     legend(beta2,'beta=0.1',' radom number')
%     title ('simulations with different beta')
%     hold off
    
end
% plot variance
figure(4)
hold on
vari_1 = var(x_1);
vari_2 = var(x_2);
plot(t,vari_1);
plot(t,vari_2);
legend('beta=0.1', beta2);
figure(5);
subplot(2,1,1);plot(t,x_1); title('beta=0.1 ');xlabel('t');ylabel('x')
subplot(2,1,2);plot(t,x_2); title( beta2 );xlabel('t');ylabel('x')

