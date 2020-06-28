% Xiao, Tianqi 3371477
% Lab2 17.11.18
% last changes 21.11.18

% This program is to plot the world cost lines as well as the parameter lines within the interval 
%       ¦«¡Ê[-180¡ã,-150¡ã,...,180¡ã], ¦µ¡Ê[0¡ã,15¡ã,...,90¡ã]
% and display the Tissot distortion ellipses together with the eigen-directions within the interval
%       ¦«¡Ê[-180¡ã,-150¡ã,...,150¡ã], ¦µ¡Ê[0¡ã,30¡ã,...,90¡ã,¦µ0]
clear all
close all
clc

% The radius of the earth is 6370 km
 R=6370;
% load data of the world coast lines and remove the point which latitude is smaller than 0
load('coast.mat');
long(long>180) = long(long>180)-360;
[long, lat] = jump2nan(long, 180,lat);
lat(lat<0)=NaN;
    % figure('color',[1,1,1]);
    % plot(long,lat)
% phi0=const, Assume phi0=45
 phi0=45;
%% plot the coast line
[x,y]=Stereographic(long,lat,phi0,R);
figure('color',[1,1,1]);hold on;
view(90,90)
plot(x,y)

axis equal
axis off

%% plot the parameter lines within the interval
% the interval of parameter lines
long_interval=-180:30:180;
lat_interval=0:15:90;
% meridians 
[long_m,lat_m]=meshgrid(long_interval ,linspace(min(lat_interval),max(lat_interval),360)); 
[x_m,y_m]=Stereographic(long_m,lat_m,phi0,R);
%check the maridians
%figure(2);
%plot(long_m,lat_m)
%figure(3)
%plot(x_m,y_m);

% parallels
[long_p,lat_p]=meshgrid(linspace(min(long_interval),max(long_interval),360),lat_interval);
[x_p,y_p]=Stereographic(long_p,lat_p,phi0,R);
%check the parallels
%figure(4);
%plot(long_p,lat_p)
%figure(5)
%plot(x_p',y_p');

%plot the parameterline
plot(x_m,y_m,'color', [0.5 0.5 0.5]);
plot(x_p',y_p','color', [0.5 0.5 0.5]);

%% Tissot ellipse

% the interval of Tissot distortion ellipses
long_interval=-180:30:150;
lat_interval=[0:30:90,phi0];

for Lamda=long_interval
    for Phi=(lat_interval)
        
        %determine the center of the ellipse
        [x0, y0] = Stereographic(Lamda, Phi, phi0 ,R);
        %fprintf('lamda=%d,delta=%d ',Lamda,Delta)
        
        %change degree into rad
        lamda=Lamda*pi/180;
        delta=(90-Phi)*pi/180;
        
       %% solve the eigenvalue problem for C and G
        %   G=[R^2*sin(delta)^2 0;0 R^2]    with delta= pi/2-phi
        %   g=[1 0;0 1];
        %   C=J'*g*J
        G=[R^2*sin(delta)^2 0;0 R^2];
        g=[1 0;0 1];
        J= [-R*(1+sin(phi0*pi/180))*tan(delta/2)*sin(lamda), R*(1+sin(phi0*pi/180))*cos(lamda)*(tan(delta/2)*tan(delta/2)/2 + 1/2);  R*(1+sin(phi0*pi/180))*tan(delta/2)*cos(lamda), R*(1+sin(phi0*pi/180))*sin(lamda)*(tan(delta/2)*tan(delta/2)/2 + 1/2)];
        C=J'*g*J; 
        % C and G are both diagonal, so the extreme distortions :
        % lamda_max^2=C11/G11=(1+sin(phi0))/(2*cos(delta/2)^2) 
        % lamda_min^2=C22/G22=(1+sin(phi0))/(2*cos(delta/2)^2)    
       
        % eigen vectors
        [F,D]=eig(C,G);
        f=J*F; 
        
        %eigen values
        lamda_max = (1+sin(phi0*pi/180))/(2*cos(delta/2)*cos(delta/2));
        lamda_min = (1+sin(phi0*pi/180))/(2*cos(delta/2)*cos(delta/2));
      
        if isnan(lamda_max)
            lamda_max=lamda_min;
        end
   
        % ang is the angle with coordinate axes
        ang = atan2(f(2,1), f(1,1));
   
        % ra and rb are the semimajor axis and the semiminor axis
        % set scale
        scale = 600;
        ra = scale * lamda_max;
        rb = scale * lamda_min;
       %% plot the Tissot ellipses
        if (abs(lamda_max-1)<1e-13 && abs(lamda_min -1)<1e-13)
        % Green circles: Tissot ellipses without distortion
            ellipse(ra, rb, ang, x0, y0, 'g'); 
        % Draw the semimajor axis and semiminor axis
            line([x0-ra*cos(ang),x0+ra*cos(ang)],[y0-ra*sin(ang),y0+ra*sin(ang)],'color','g');
            line([x0-rb*sin(pi-ang),x0+rb*sin(pi-ang)],[y0-rb*cos(pi-ang),y0+rb*cos(pi-ang)],'color','g');
 
        else
        % Red circles: Tissot ellipses with distortion in different direction
            ellipse(ra, rb, ang, x0, y0, 'r'); 
        % Draw the semimajor axis and semiminor axis
            line([x0-ra*cos(ang),x0+ra*cos(ang)],[y0-ra*sin(ang),y0+ra*sin(ang)],'color','r');
            line([x0-rb*sin(pi-ang),x0+rb*sin(pi-ang)],[y0-rb*cos(pi-ang),y0+rb*cos(pi-ang)],'color','r');
        end
    end
end
title('"Stereographic Projection"')
print('-depsc2', 'figure.eps');
print('-dpng', '-r400', 'figure.png')

hold off;
