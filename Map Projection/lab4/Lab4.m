% Lab4 
% Xiao Tianqi 3371477
% Wang Zhenqiao 3371590
% 25.01.19
clear
clc

%% ellipsoid coordinates of Hoellbuck(LH,BH) , Kapfwald(LK,BK) and Local origin P0

% point Kapfwald(LK,BK)
BK_dms=[48 01 18.2679];
LK_dms=[8 33 07.9143];

% point Hoellbuck(LH,BH)
LH_dms=[10 16 57.7583];
BH_dms=[48 43 37.8185];% Reference maridian L0 = 9 degree

% Local origin P0 (L0,B0) consist in the integer part of points of interest, 
% B0 = 48 degree
L0_dms=[9 0 0];
B0_dms=[48 0 0];

% dms to rad 
L0_rad=deg2rad(dms2degrees(L0_dms));
B0_rad=deg2rad(dms2degrees(B0_dms));
LK_rad=deg2rad(dms2degrees(LK_dms));
BK_rad=deg2rad(dms2degrees(BK_dms));
LH_rad=deg2rad(dms2degrees(LH_dms));
BH_rad=deg2rad(dms2degrees(BH_dms));

% GRS80 Ellipsoid
A= 6378137; % [m]
E_2 =0.00669438002;

% Bassel Ellipsoid
A_B= 6377397.155; % [m]
E_2B=0.0066743722;

% b= Bp - B0; l= Lp - L0
bK = BK_rad-B0_rad; % [rad]
bH = BH_rad-B0_rad; % [rad]
lK = LK_rad-L0_rad; % [rad]
lH = LH_rad-L0_rad; % [rad]
fprintf('Transform from ellipsoidal coordinates to different Cartesian coordinates\r\n')
%% the Gauss-Krueger coordinates 
m0_GK=1;
fprintf('Gauss-Krueger Projection:\r\n')
fprintf('Point Hoellbuck:\r\n')
% conformal coordinates of point Hoellbuck (LH,BH) => x(l,b), y(l,b)  
GK_Hoe=Cal_xy(m0_GK,lH,bH,B0_rad);
fprintf('The conformal coordinates (x,y):\r\n   (%.4f m ,%.4f m)\r\n',GK_Hoe(1),GK_Hoe(2));

% Northing X(l,B) = X0 + x(l,b), with X0=lat2mer(A,E_2,B0)
GK_Hoe_X= GK_Hoe(1)+ m0_GK*lat2mer(A,E_2,B0_rad);
% False Easting Y = y(l.b) + (L0/3)*10e6 + 500000
GK_Hoe_Y= GK_Hoe(2)+ dms2degrees(L0_dms)/3*10e6 +500000;
fprintf('The Gauss-Krueger coordinats:\r\n   Northing H = %.4f m , False Easting R = %.4f m \r\n',GK_Hoe_X,GK_Hoe_Y);
% meridian convergence and distortion
GK_Hoe_c= rad2deg(GK_Hoe(3)); % degree
GK_Hoe_cdms=degrees2dms(GK_Hoe_c); % dms
GK_Hoe_lamda= GK_Hoe(4);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',GK_Hoe_c,GK_Hoe_cdms,GK_Hoe_lamda);

% point Kapfwald(LK,BK)
fprintf('\r\nPoint Kapfwald:\r\n')
GK_Kap=Cal_xy(m0_GK,lK,bK,B0_rad);
fprintf('The conformal coordinates (x,y):\r\n   (%.4f m ,%.4f m)\r\n',GK_Kap(1),GK_Kap(2));
% Northing and False Easting
GK_Kap_X= GK_Kap(1)+ m0_GK*lat2mer(A,E_2,B0_rad);
GK_Kap_Y= GK_Kap(2)+  dms2degrees(L0_dms)/3*10e6 +500000;
fprintf('The Gauss-Krueger coordinats:\r\n   Northing H = %.4f m , False Easting R = %.4f m \r\n',GK_Kap_X,GK_Kap_Y);
% meridian convergence and distortion
GK_Kap_c=rad2deg(GK_Kap(3)); %degree
GK_Kap_cdms=degrees2dms(GK_Kap_c); %dms
GK_Kap_lamda= GK_Kap(4);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',GK_Kap_c,GK_Kap_cdms,GK_Kap_lamda);

%% the UTM coordinates 
% check
% b=7.78115826e-3
% l=2.995461569e-2
%  (l,b) => (lK,bK) and (lH ,bH) 
m0_UTM=0.9996;
fprintf('\r\nUTM Projection:\r\n')
fprintf('Point Hoellbuck:\r\n')
% point Hoellbuck(LH,BH)
UTM_Hoe=Cal_xy(m0_UTM,lH,bH,B0_rad);
fprintf('The conformal coordinates (x,y):\r\n   (%.4f m ,%.4f m)\r\n',UTM_Hoe(1),UTM_Hoe(2));
% False Northing and False Easting 
UTM_Hoe_X= UTM_Hoe(1)+ m0_UTM*lat2mer(A,E_2,B0_rad);
UTM_Hoe_Y= UTM_Hoe(2)+500000;
fprintf('The UTM coordinats:\r\n   False Northing N = %.4f m , False Easting E = %.4f m \r\n',UTM_Hoe_X,UTM_Hoe_Y);
% meridian convergence and distortion
UTM_Hoe_c=rad2deg(UTM_Hoe(3));
UTM_Hoe_cdms=degrees2dms(UTM_Hoe_c);
UTM_Hoe_lamda= UTM_Hoe(4);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',UTM_Hoe_c,UTM_Hoe_cdms,UTM_Hoe_lamda);

% point Kapfwald(LK,BK)
fprintf('\r\nPoint Kapfwald:\r\n')
UTM_Kap=Cal_xy(m0_UTM,lK,bK,B0_rad);
fprintf('The conformal coordinates (x,y):\r\n   (%.4f m ,%.4f m)\r\n',UTM_Kap(1),UTM_Kap(2));
% False Northing and False Easting 
UTM_Kap_X= UTM_Kap(1)+ m0_UTM*lat2mer(A,E_2,B0_rad);
UTM_Kap_Y= UTM_Kap(2)+500000;
fprintf('The UTM coordinats:\r\n   False Northing N = %.4f m , False Easting E = %.4f m \r\n',UTM_Kap_X,UTM_Kap_Y);
% meridian convergence and distortion
UTM_Kap_c= rad2deg(UTM_Kap(3));
UTM_Kap_cdms=degrees2dms(UTM_Kap_c);
UTM_Kap_lamda= UTM_Kap(4);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',UTM_Kap_c,UTM_Kap_cdms,UTM_Kap_lamda);

%% Transback
fprintf('\r\nTransform from Cartesian coordinates to ellipsoidal coordinates in respect of Bessel Ellipsoid\r\n')
X0_2=5400000;
% UTM coordinates
fprintf('UTM coordinates transback to Bessel ellipsoid\r\n')
% check
% x=-32723.2852
% y=126917.0456
B0_UTM=mer2lat(A_B,E_2B,X0_2/m0_UTM);

% point Hoellbuck(xH,yH)
fprintf('Point Hoellbuck:\r\n')
% compute l(x,y) and b(x,y) c and lamda
UTM_Hoe_x= UTM_Hoe_X - X0_2;
UTM_Hoe_y= UTM_Hoe_Y - 500000;
UTM_Ell_points_H = Cal_lb(m0_UTM,UTM_Hoe_x,UTM_Hoe_y,B0_UTM);
% compute (Lp,Bp)
UTM_Hoe_B_deg=rad2deg(B0_UTM+UTM_Ell_points_H(1));
UTM_Hoe_B_dms=degrees2dms(UTM_Hoe_B_deg);
UTM_Hoe_L_deg=rad2deg(L0_rad+UTM_Ell_points_H(2));
UTM_Hoe_L_dms=degrees2dms(UTM_Hoe_L_deg);
fprintf('Bessel-ellipsoidal coordinats:\r\nIn decimal degrees: B = %.4f¡ã, L = %.4f¡ã\r\nIn dms: B = %d¡ã%d¡¯%.4f" , L = %d¡ã%d¡¯%.4f" \r\n',UTM_Hoe_B_deg,UTM_Hoe_L_deg,UTM_Hoe_B_dms,UTM_Hoe_L_dms);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',rad2deg(UTM_Ell_points_H(3)),degrees2dms(rad2deg(UTM_Ell_points_H(3))),UTM_Ell_points_H(4));

% point Kapfwald (xK,yK)
fprintf('\r\nPoint Kapfwald:\r\n')
% compute l(x,y) and b(x,y) c and lamda
UTM_Kap_x= UTM_Kap_X - X0_2;
UTM_Kap_y= UTM_Kap_Y - 500000;
UTM_Ell_points_K = Cal_lb(m0_UTM,UTM_Kap_x,UTM_Kap_y,B0_UTM);
% compute (Lp,Bp)
UTM_Kap_B_deg=rad2deg(B0_UTM+UTM_Ell_points_K(1));
UTM_Kap_B_dms=degrees2dms(UTM_Kap_B_deg);
UTM_Kap_L_deg=rad2deg(L0_rad+UTM_Ell_points_K(2));
UTM_Kap_L_dms=degrees2dms(UTM_Kap_L_deg);
fprintf('Bessel-ellipsoidal coordinats:\r\nIn decimal degrees: B = %.4f¡ã, L = %.4f¡ã\r\nIn dms: B = %d¡ã%d¡¯%.4f" , L = %d¡ã%d¡¯%.4f" \r\n',UTM_Kap_B_deg,UTM_Kap_L_deg,UTM_Kap_B_dms,UTM_Kap_L_dms);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',rad2deg(UTM_Ell_points_K(3)),degrees2dms(rad2deg(UTM_Ell_points_K(3))),UTM_Ell_points_K(4));

% Gauss-Krueger coordinates
fprintf('\r\nGauss-Krueger coordinates transback to Bessel ellipsoid\r\n')
B0_GK=mer2lat(A_B,E_2B,X0_2/m0_GK);

% point Hoellbuck(xH,yH)
fprintf('Point Hoellbuck:\r\n')
% compute l(x,y) and b(x,y) c and lamda
GK_Hoe_x= GK_Hoe_X - X0_2;
GK_Hoe_y= GK_Hoe_Y - 500000- dms2degrees(L0_dms)/3*10e6;
GK_Ell_points_H = Cal_lb(m0_GK,GK_Hoe_x,GK_Hoe_y,B0_GK);

% compute Lp and Bp
GK_Hoe_B_deg=rad2deg(B0_GK+GK_Ell_points_H(1));
GK_Hoe_B_dms=degrees2dms(GK_Hoe_B_deg);
GK_Hoe_L_deg=rad2deg(L0_rad+GK_Ell_points_H(2));
GK_Hoe_L_dms=degrees2dms(GK_Hoe_L_deg);
fprintf('Bessel-ellipsoidal coordinats:\r\nIn decimal degrees: B = %.4f , L = %.4f¡ã\r\nIn dms: B = %d¡ã%d¡¯%.4f" , L = %d¡ã%d¡¯%.4f" \r\n',GK_Hoe_B_deg,GK_Hoe_L_deg,GK_Hoe_B_dms,GK_Hoe_L_dms);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',rad2deg(GK_Ell_points_H(3)),degrees2dms(rad2deg(GK_Ell_points_H(3))),GK_Ell_points_H(4));

% point Kapfwald (xK,yK)
fprintf('\r\nPoint Kapfwald:\r\n')
% compute l(x,y) and b(x,y) c and lamda
GK_Kap_x= GK_Kap_X - X0_2;
GK_Kap_y= GK_Kap_Y - 500000- dms2degrees(L0_dms)/3*10e6;
GK_Ell_points_K = Cal_lb(m0_GK,GK_Kap_x,GK_Kap_y,B0_GK);
% compute Lp and Bp
GK_Kap_B_deg=rad2deg(B0_GK+GK_Ell_points_K(1));
GK_Kap_B_dms=degrees2dms(GK_Kap_B_deg);
GK_Kap_L_deg=rad2deg(L0_rad+GK_Ell_points_K(2));
GK_Kap_L_dms=degrees2dms(GK_Kap_L_deg);
fprintf('Bessel-ellipsoidal coordinats:\r\nIn decimal degrees: B = %.4f , L = %.4f¡ã\r\nIn dms: B = %d¡ã%d¡¯%.4f" , L = %d¡ã%d¡¯%.4f" \r\n',GK_Kap_B_deg,GK_Kap_L_deg,GK_Kap_B_dms,GK_Kap_L_dms);
fprintf('The meridian convergence:  c[degrees] = %.4f¡ã  c[dms] = %d¡ã%d¡¯%.4f" \r\nDistortion:   Lamda = %.4f\r\n',rad2deg(GK_Ell_points_K(3)),degrees2dms(rad2deg(GK_Ell_points_K(3))),GK_Ell_points_K(4));

%% comparison of different reference ellipsoid
% Gauss-Kruger Projection
% delta = Bessel coordinates - GRS80 coordinates
delta_GK =[GK_Kap_B_deg-rad2deg(BK_rad) GK_Hoe_B_deg-rad2deg(BH_rad); GK_Kap_L_deg-rad2deg(LK_rad) GK_Hoe_L_deg-rad2deg(LH_rad)];
fprintf('\r\nInfluence due to the changed geometry of Gauss-Krueger coordinates in second: \r\nKapfwald:  delta B = %.4f" delta L = %.4f" \r\nHoellbuck: delta B = %.4f" delta L = %.4f"\r\n ',3600*delta_GK);
% UTM Projection
% delta = Bessel coordinates - GRS80 coordinates
delta_UTM=[UTM_Kap_B_deg-rad2deg(BK_rad) UTM_Hoe_B_deg-rad2deg(BH_rad); UTM_Kap_L_deg-rad2deg(LK_rad) UTM_Hoe_L_deg-rad2deg(LH_rad)];
fprintf('\r\nInfluence due to the changed geometry of UTM coordinates in second:  \r\nKapfwald:  delta B = %.4f" delta L = %.4f" \r\nHoellbuck: delta B = %.4f" delta L = %.4f"\r\n',3600*delta_UTM);

% approximate estimation
% delta S = G(B_bessel)-G(B_grs)
fprintf('\r\nApproximative estimation of the "influence" in meter: \r\n')
% Gauss-Krueger Projection
deltaS_GK  =[lat2mer(A_B,E_2B,B0_GK+GK_Ell_points_K(1)) - lat2mer(A,E_2,BK_rad);lat2mer(A_B,E_2B,B0_GK+GK_Ell_points_H(1)) - lat2mer(A,E_2,BH_rad)];
fprintf('Transback of Gauss-Krueger coordinates:    Kapfwald: delta S = %.4f m 	 Hoellbuck:delta S = %.4f m\r\n',deltaS_GK);
% UTM Projection
deltaS_UTM =[lat2mer(A_B,E_2B,B0_UTM+UTM_Ell_points_K(1))-lat2mer(A,E_2,BK_rad);lat2mer(A_B,E_2B,B0_UTM+UTM_Ell_points_H(1)) - lat2mer(A,E_2,BH_rad)];
fprintf('Transback of UTM coordinates:              Kapfwald: delta S = %.4f m 	 Hoellbuck:delta S = %.4f m\r\n',deltaS_UTM);
