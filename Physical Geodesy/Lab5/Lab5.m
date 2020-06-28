% Physical Geodesy Assignment 5
% Xiao, Tianqi
% 3371477
% 11.07.2019

% clear all

theta=0:180;
%% Task 1
l=10;
m=0:10;
syms t
% Normalized Legnendre function Rodrigues-Ferrers
pl=1/(2^l*factorial(l))*diff((t^2-1)^l,t,l);    % Pl(t)

for i=1:length(m)
if m(i)==0
    N=sqrt(2*l+1);
else
    N=sqrt(2*(2*l+1)*(factorial(l-m(i))/factorial(l+m(i))));
end
Plm_RF(i,:)=double(subs((1-t^2).^(m(i)/2).*(diff(pl,t,m(i)))*N,t,cos(deg2rad(theta)))); % Plm(t)
labely_RF{i}=['m=' num2str(m(i))];
end


% Normalized Legnendre function Recursive Formulas

% initials
P(1,1)=t.^0; % P0,0
P(2,1)=sqrt(3)*t; % P1,0
P(2,2)=sqrt(3*(1-t.^2));% P1,1
for i=2:l % l 
    P(i+1,i+1)=sqrt((2*i+1)/(2*i)*(1-t.^2))*P(i,i); % Pl,l
    for j=0:i-1 % m
        P(i+1,j+1)=sqrt((2*i+1)/((i+j)*(i-j)))*(sqrt(2*i-1)*t*P(i,j+1)-sqrt((i-1+j)*(i-1-j)/(2*i-3))*P(i-1,j+1)); % Pl,m
     if i==l && ismember(j,m)==1
        Plm_Re(find(m==j),:)=double(subs(P(i+1,j+1),t,cos(deg2rad(theta)))); % Plm(t)
        labely_Re{find(m==j)}=['m=' num2str(j)];
     end
    end
    if i==l
    Plm_Re(find (m==i),:)=double(subs(P(i+1,i+1),t,cos(deg2rad(theta)))); % Plm(t)
    labely_Re{find (m==i)}=['m=' num2str(i)];
    end
end
% Plot
figure (1)
subplot(1,2,1)
axis normal
grid on
hold on
plot(theta,m'.*4+Plm_RF)
yticks(m*4)
yticklabels(labely_RF)
xticks(0:30:180);
xlabel('\theta [deg]')
title({'Legendre Funcion of degree 10 ','Using Rodrigues-Ferrers'})

subplot(1,2,2)
grid on
hold on
plot(theta,m'.*4+Plm_Re)
yticks(m*4)
yticklabels(labely_Re)
xticks(0:30:180);
xlabel('\theta [deg]')
title({'Legendre Funcion of degree 10 ','Using  Recursive Formulas'})

% Spherical Harmony
lambda=-180:0.5:180;    %longtitude

% l = 10 Calculate Y(theta,lamda) for different m
% figure in plane
figure 
for i= 1:length(m)
Y(:,:,i)=Plm_RF(i,:)'.*cos(m(i)*deg2rad(lambda));
subplot(3,4,i)
hold on
imagesc(theta,lambda,Y(:,:,i));
yticks(-180:60:180);
xticks(0:30:180);
colorbar;
title (['Y(\theta,\lambda) with l =',num2str(l),'m = ',num2str(m(i))]);
end
colormap(jet);

% Figure in sphere
r=5;
[lat,long]=meshgrid(deg2rad(90-theta),deg2rad(lambda));
[x,y,z] = sph2cart(long,lat,r);
figure
for i=1:length(m)
subplot(3,4,i)
hold on
axis equal
surf(x,y,z,Y(:,:,i)','EdgeColor', 'none');
colormap(jet);
colorbar;
title (['Y(\theta,\lambda) with l=10 , m= ',num2str(i-1)])
view(3)
end
%%  Task 2
t_p=90;
t_q=0:90;
% psi_PQ
psi_pq=sin(deg2rad(t_p)).*sin(deg2rad(t_q))+cos(deg2rad(t_p)).*cos(deg2rad(t_q)); 
% LHS
P_left(1,:)=psi_pq.^0; % P0,0
P_left(2,:)=psi_pq; % P1,0

% RHS
P_tp(1,1,:)=cos(deg2rad(t_p)).^0; % P0,0
P_tp(2,1,:)=sqrt(3)*cos(deg2rad(t_p)); % P1,0
P_tp(2,2,:)=sqrt(3*(1-cos(deg2rad(t_p)).^2));% P1,1

P_tq(1,1,:)=cos(deg2rad(t_q)).^0; % P0,0
P_tq(2,1,:)=sqrt(3)*cos(deg2rad(t_q)); % P1,0
P_tq(2,2,:)=sqrt(3*(1-cos(deg2rad(t_q)).^2));% P1,1

P_right(1,:)=P_tp(1,1,:)*P_tq(1,1,:);
P_right(2,:)=(sum(P_tp(2,:,:).*P_tq(2,:,:)))/3;

for l=2:100             % l=0:100
    % get LHS value
    P_left(l+1,:)=((2*l-1)*psi_pq.*P_left(l,:)-(l-1)*P_left(l-1,:))/l;
    for m=0:l % m
        % get p_lm(cosp) and p_lm(cosq)
        if m==l         % Pl,l
            P_tp(l+1,m+1,:)=sqrt((2*l+1)/(2*l))*(P_tp(2,2,:)/sqrt(3)).*P_tp(l,m,:);
            P_tq(l+1,m+1,:)=sqrt((2*l+1)/(2*l))*(P_tq(2,2,:)/sqrt(3)).*P_tq(l,m,:);
        else
            P_tp(l+1,m+1,:)=(sqrt((2* l-1))*(P_tp(2,1,:)/sqrt(3)).*P_tp(l,m+1,:)-sqrt((l-1+m)*(l-1-m)/(2*l-3))*P_tp(l-1,m+1,:))*sqrt((2*l+1)/(l+m)/(l-m));
            P_tq(l+1,m+1,:)=(sqrt((2* l-1))*(P_tq(2,1,:)/sqrt(3)).*P_tq(l,m+1,:)-sqrt((l-1+m)*(l-1-m)/(2*l-3)).*P_tq(l-1,m+1,:))*sqrt((2*l+1)/(l+m)/(l-m));
        end
    end
    % RHS
    P_right(l+1,:)=(sum(P_tp(l+1,:,:).*P_tq(l+1,:,:)))/(2*l+1);
end

dp=P_right-P_left;          % difference
% Plot differnce
figure
hold on;
    for l = 0:1:100
    plot(rad2deg(acos(psi_pq)), dp(l+1,:));
    end
title('Difference for different \psi with l = 0 ~100');
xlabel('Spherical distance \psi_P_Q [deg]');
ylabel('Difference');
hold off;

figure
hold on
for i = 0:1:90
    plot(0:100, dp(:,i+1));
end
title('Difference for different l (\psi = 0 - 90 )');
xlabel('Degree l');
ylabel('Difference');
hold off;

figure
hold on
surf(90:-1:0,0:100,dp,'EdgeColor','none');
colormap(jet);
colorbar;
view(3)
title('Difference for different l and \psi )');
xlabel('\psi');
ylabel('l');
zlabel('Difference');

%% Task 3 
% theta_p=theta_q=theta

% RHS
P_t(1,1,:)=cos(deg2rad(theta)).^0; % P0,0
P_t(2,1,:)=sqrt(3)*cos(deg2rad(theta)); % P1,0
P_t(2,2,:)=sqrt(3*(1-cos(deg2rad(theta)).^2));% P1,1

P_RHS(1,:)=P_t(1,1,:).*P_t(1,1,:);
P_RHS(2,:)=(P_t(2,1,:).*P_t(2,1,:)+P_t(2,2,:).*P_t(2,2,:))/3;

for l=2:100             % l=0:100
    for m=0:l % m
        % get p_lm(cosp) and p_lm(cosq)
        if m==l         % Pl,l
            P_t(l+1,m+1,:)=sqrt((2*l+1)/(2*l))*(P_t(2,2,:)/sqrt(3)).*P_t(l,m,:);
        else
            P_t(l+1,m+1,:)=(sqrt(2*l-1)*(P_t(2,1,:)/sqrt(3)).*P_t(l,m+1,:)-sqrt((l-1+m)*(l-1-m)/(2*l-3))*P_t(l-1,m+1,:))*sqrt((2*l+1)/(l+m)/(l-m));
        end
    end
    P_RHS(l+1,:)=(sum(P_t(l+1,:,:).*P_t(l+1,:,:)))/(2*l+1);
end

% Plot result
figure  
hold on;
for i = 0:180
    plot(0:100, P_RHS(:,i+1));
end
title('The right hand side of the equation for different degree');
xlabel('Degree l');
ylabel('Polynomial value');
hold off;

figure  
hold on;
for i = 0:100
    plot(0:180, P_RHS(i+1,:));
end
title('The right hand side of the equation for different \theta');
xlabel('\theta [deg]');
ylabel('Polynomial value');
hold off;

figure
hold on
surf(180:-1:0,0:100,P_RHS,'EdgeColor','none');
colormap(jet);
colorbar;
view(3)
title('The right hand side of the equation for different l and \psi )');
xlabel('\theta');
ylabel('l');
zlabel('Polynomial value');

%% Task 4

k = 11;
R = 6378136.3;              % Radius [m]
GM = 3.986004415*1e14;      % [m^3s-2]
omega = 7.292115*1e-5;      % [s-5]


r = 6379245.458;            % [m]
Lambda_P=deg2rad(10+k);     % [rad]
Theta_P = deg2rad(42+k);    % [rad]
data = load('EGM96.txt');
l_P = data(:,1);
l_max = max(l_P);
m_P = data(:,2);
clm = data(:,3);
slm = data(:,4);

% right
Plm_P(1,1)=cos(Theta_P).^0; % P0,0
Plm_P(2,1)=sqrt(3)*cos(Theta_P); % P1,0
Plm_P(2,2)=sqrt(3*(1-cos(Theta_P).^2));% P1,1
Sum_l=0;
for l=0:max(l_P) % l=0:100
    Sum_m=0;
    for m=0:l % m
        Clm_P(l+1,m+1)=clm(intersect(find(l_P==l ), find(m_P==m)));
        Slm_P(l+1,m+1)=slm(intersect(find(l_P==l ), find(m_P==m)));
        if m==l && l>1      % Pl,l
            Plm_P(l+1,m+1)=sqrt((2*l+1)/(2*l))*(Plm_P(2,2)/sqrt(3)).*Plm_P(l,m);
        else
            if l>1
                Plm_P(l+1,m+1)=(sqrt(2*l-1)*(Plm_P(2,1)/sqrt(3)).*Plm_P(l,m+1)-sqrt((l-1+m)*(l-1-m)/(2*l-3))*Plm_P(l-1,m+1))*sqrt((2*l+1)/(l+m)/(l-m));
            end
        end
        Sum_m=Sum_m+ Plm_P(l+1,m+1)*(Clm_P(l+1,m+1)*cos(m*Lambda_P)+Slm_P(l+1,m+1)*sin(m*Lambda_P));
    end
    Sum_l=Sum_l+Sum_m*(R/r)^(l+1);
end
V=GM/R*Sum_l
Vc = omega^2*r^2*sin(Theta_P)^2/2
W = Vc + V