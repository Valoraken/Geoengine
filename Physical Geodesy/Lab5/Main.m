% Physical Geodesy Lab 5
% Name : Yi Hong/Shu Suo
% Matriculation number: 3294211/3295472
%
clc

%% Task 1
% Recursive method
theta = 0:1:180;
t = cos(deg2rad(theta));
l = 10;
Plm_rec = NormLegendreFunc(l, t);

% Rodrigues-Ferrers method
syms x;
Plm_RF= zeros(11,11,181);
Plm_RF(1,1,:)=ones(1,181);
for i = 1:l
    y=(x.^2-1).^i;
    dy = diff(y,i)/(2^i*factorial(i));
    Plm_RF(i+1,1,:) = sqrt(2*i+1)*subs(dy,x,t);
    for j = 1:i
        ddy = (1-x^2)^(j/2).*diff(dy,j);
        Plm_RF(i+1,j+1,:) = sqrt(2*(2*i+1)*factorial(i-j)/factorial(i+j))*subs(ddy,x,t); 
    end
end
d_t1 = Plm_rec-Plm_RF;
d1 = zeros(1,181);
figure(1);
for l = 0:1:5
    d1(l+1,:) = reshape(d_t1(11,l+1,:),[],1);
    plot(theta, d1(l+1,:));
    hold on;
end
xlabel('Colatitude (?)');
ylabel('Difference');
title('The Difference between Recursion and Rodrigues-Ferrers (Degree:10,Order:0-10)');
hold off;

figure(2);
P10_1 = zeros(6, 181);
for l = 0:1:5
    P10_1(l+1,:) = reshape(Plm_RF(11,l+1,:),[],1);
    plot(90-theta, P10_1(l+1,:));
    hold on;
end
legend('P_{10,0}','P_{10,1}','P_{10,2}','P_{10,3}','P_{10,4}','P_{10,5}');
xlabel('Latitude \phi (?)');
ylabel('Normalized Legendre function');
title('Fully normalized Legendre Function of Degree 10(order 0-5)');
axis([-90 90 -5 5]);
hold off;

figure(3);
P10_2 = zeros(5, 181);
for l = 6:1:10
    P10_2(l-5,:) = reshape(Plm_RF(11,l+1,:),[],1);
    plot(90-theta, P10_2(l-5,:));
    hold on;
end
legend('P_{10,6}','P_{10,7}','P_{10,8}','P_{10,9}','P_{10,10}');
title('Fully normalized Legendre Function of Degree 10(order 6-10)');
xlabel('Latitude \phi (?)');
ylabel('Normalized Legendre function');
axis([-90 90 -3 3]);
hold off;
% 
lambda = 0:1:360; 
Y = zeros(181,361,11);
P10 = [P10_1;P10_2];
for m = 0:1:10
    Q = cos(m*deg2rad(lambda));
    [q, p] = meshgrid(Q, P10(m+1,:));
    Y(:,:,m+1) = q.*p;
end
figure(4);
[lam,phi] = meshgrid(deg2rad(lambda),deg2rad(90-theta));
R = ones(size(lam));
[x,y,z] = sph2cart(lam, phi, R);
surf(x,y,z,Y(:,:,6),'EdgeColor', 'none');
axis equal;
colormap('jet');
% planar visualization
% surf(lambda,theta,Y(:,:,6), 'EdgeColor', 'none');
% view(2);
% axis([lambda(1) lambda(end) theta(1) theta(end)]);

%% Task 2
% Legendre polynomial from spherical distance 
theta_p = 90; 
t_p = cos(deg2rad(theta_p)); 
Phi_p = deg2rad(90-theta_p); 
theta_q = 0:1:90; 
t_q = cos(deg2rad(theta_q)); 
Phi_q = deg2rad(90-theta_q); 
t_pq = sin(Phi_p)*sin(Phi_q)+cos(Phi_p)*cos(Phi_q); 
PHI = acos(t_pq)*180/pi;
P1 = LegendreFunc(100,t_pq);
P_pq1 = reshape(P1(:,1,:),101,91);

% Legendre polynomial from individual contribution 
P_p = NormLegendreFunc(100, t_p); 
P_q = NormLegendreFunc(100, t_q); 
P_pq2 = zeros(91,101); 
c = 0; 
for l = 0:1:100
    temp = zeros(1,1,91); 
    for m = 0:1:l
        temp = temp + P_p(l+1,m+1,:).*P_q(l+1,m+1,:);
    end
    temp = temp/(2*l+1); P_pq2(:,l+1) = temp(1,1,:); c = c+1;
end
d_t2 = P_pq2'- P_pq1;

figure(5);
for l = 0:1:100
    plot(PHI, d_t2(l+1,:));
    hold on;
end
title('Difference for different spherical distance (degree:0-100)');
xlabel('Spherical distance \Phi (?)');
ylabel('Difference');
hold off;

figure(6);
l_t2 = 0:1:100;
for sd = 0:1:90
    plot(l_t2, d_t2(:, sd+1));
    hold on;
end
title('Difference for different degree (spherical distance: 0?-90?)');
xlabel('Degree l');
ylabel('Difference');
hold off;

%% Task 3
P2 = NormLegendreFunc(100, t);
P_t3 = zeros(181,101);
c = 0;
for l = 0:1:100
    temp = zeros(1,1,181);
    for m = 0:1:l 
        temp = temp + P2(l+1,m+1,:).^2;  
    end
    temp = temp/(2*l+1);
    P_t3(:,l+1) = reshape(temp(1,1,:),181,1);
    c = c+1;
end

figure(7);
l_t3 = 0:1:100;
for theta = 0:1:180
    plot(l_t3, P_t3(theta+1,:));
    hold on;
end
title('The right hand side of the equation for different degree');
xlabel('Degree l');
ylabel('Polynomial value');
hold off;

%% Task 4
k = 11;
R = 6378136.3;
GM = 3.986004415*10^(14);
omega = 7.292115*10^(-5);
r = 6379245.458;

lambda_P = deg2rad(10+k);
theta_P = deg2rad(42+k);
Data = load('EGM96.txt');
l_t4 = Data(:,1);
l_max = max(l_t4);
m_t4 = Data(:,2);
clm = Data(:,3);
slm = Data(:,4);

t_P = cos(theta_P);
P_P = NormLegendreFunc(l_max, t_P);
P_P = reshape(P_P',[],1);
P_P(P_P==0)=[];
temp = P_P.*(clm.*cos(m_t4*lambda_P)+slm.*sin(m_t4*lambda_P));
idx = 1:1:l_max+1;
temp1 = mat2cell(temp,idx);
V_P = 0;
VPa=[]
for idx = 1:1:l_max+1
    V_P = V_P+(R/r)^idx*sum(cell2mat(temp1(idx)));
    VPa=[VPa sum(cell2mat(temp1(idx)))];
end
V_P = GM/R*V_P;
Vc_P = omega^2*r^2*sin(theta_P)^2/2;
W_P = Vc_P + V_P;

%% Task 5 
r1 = R:100:r;
V_P1 = zeros(1,size(r1,2));
Vc_P1 = zeros(1,size(r1,2));
for i = 1:length(r1)
    Vc_P1(i) = omega^2*r1(i)^2*sin(theta_P)^2/2;
    for idx = 1:1:l_max+1
        V_P1(i) = V_P1(i)+(R/r1(i))^idx*sum(cell2mat(temp1(idx)));
    end
end
V_P1 = GM/R*V_P1;
W_P1 = Vc_P1 + V_P1;
figure(8);
plot(r1, Vc_P1);
title('The centrifugal potential varies with different height');
xlabel('Distace to Earth center (m)');
ylabel('Ccentrifugal potential (m^2/s^2)');
figure(9);
plot(r1, V_P1);
hold on;
plot(r1, W_P1);
title('The gravity and gravitational potential varies with different height');
xlabel('Distace to Earth center (m)');
ylabel('Potential (m^2/s^2)');
legend('Gravitational potential','Gravity potential');
hold off;