run main.m
%% Define Parameters
GM = 3.98600415*1e14;
t = 561600
OmegaDOT_e = 7.2921151467*1e-5  % WGS84 Earth rotation rate.
%% Compute for Satellight 14
X_14=getGPS(SV, 14, t);
%% Compute for Satellite 18
X_18=getGPS(SV, 18, t);
Dis14_18=sqrt(sum((X_18-X_14).^2))


%% Compute for Satellite 10
t=564600

[X_10 V_10]=getGPS(SV, 10, t)

%% Compute for Satellite 10 (t=t0-1)
t1 = t-1
X_10_1=getGPS(SV, 10, t1)

%% Compute for Satellite 10_2_2
t2 = t+1
X_10_2=getGPS(SV, 10, t2)
v_10=(X_10_2-X_10_1)/2

diff_v_10=V_10-v_10