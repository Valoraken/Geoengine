function [X,VX] = getGPS(SV,n,t)
% n satellite ID
% t TOW 
GM = 3.98600415*1e14;
OmegaDOT_e = 7.2921151467*1e-5;  % WGS84 Earth rotation rate.
%% Compute for Satellight n
a0= (SV(n).navData.sqrtA)^2;
M0= SV(n).navData.M0;
n0= sqrt(GM/a0^3);
t0= SV(n).navData.TOE;
dn= SV(n).navData.DeltaN;
M= M0+(n0+dn)*(t-t0);
e=SV(n).navData.e;
%Normalize mean anomaly
while max(M>= 2*pi)
    M(M>= 2*pi) = M(M>= 2*pi)-2*pi; 
end
while max(M< 0)    
    M(M< 0)     = M(M< 0)    +2*pi; 
end

% iteration to get eccentric anomaly
E0=M; % E0
E= M+e.*sin(E0);
Ei=E0;
while max(abs(Ei-E))>1e-14
    Ei= E;
    E= M+e.*sin(Ei);
end
niu= 2*atan2(sqrt((1+e)./(1-e)),1./tan(E./2));
omega=SV(n).navData.omega;
Cuc=SV(n).navData.Cuc;
Cus=SV(n).navData.Cus;
u=omega+niu+Cuc*cos(2*(omega+niu))+Cus*sin(2*(omega+niu)); % u(t)
Crc=SV(n).navData.Crc;
Crs=SV(n).navData.Crs;
r=a0*(1-e*cos(E))+Crc*cos(2*(omega+niu))+Crs*sin(2*(omega+niu));
i0=SV(n).navData.i0;
IDOT=SV(n).navData.IDOT;
Cic=SV(n).navData.Cic;
Cis=SV(n).navData.Cis;
i=i0+IDOT*(t-t0)+Cic*cos(2*(omega+niu))+Cis*sin(2*(omega+niu)); % i(t)
Omega0=SV(n).navData.OMEGA0;
OmegaDOT=SV(n).navData.OMEGA_DOT;
OMEGA=Omega0+(OmegaDOT-OmegaDOT_e)*(t-t0)-OmegaDOT_e*t0;          %  lambda(t)

X=angle2dcm( -u, -i, -OMEGA, 'ZXZ' )*[r;0;0];

%% Velocity
MDOT=n0+dn;
EDOT= MDOT/(1-e*cos(E));
%???
% niuDOT2= sin(E)*EDOT*(1.0+e*cos(niu))/(sin(niu)*(1-e*cos(E)));
niuDOT= sqrt((1+e)./(1-e))*((cos(niu/2)^2)/(cos(E/2)^2))*EDOT;
OMEGADOT=OmegaDOT-OmegaDOT_e;
niuDOT3= a0*sqrt(1-e^2)*(n+n0)/r^2;

rDOT= a0*e*sin(E)*EDOT-2*((Crc*cos(2*(niu+omega))+Crs*sin(2*(niu+omega)))*niuDOT);
iDOT=IDOT+2*(Cic*cos(2*(niu+omega))+Cis*sin(2*(niu+omega)))*niuDOT;
uDOT=niuDOT+2*(Cus*cos(2*(niu+omega))-Cuc*sin(2*(niu+omega)))*niuDOT;
% uDOT2=niuDOT2+2*(Cus*cos(2*u)-Cuc*sin(2*u))*niuDOT2;
% rDOT2=a0^2*e*sin(E)*(n0+dn)/r+2*(Crc*cos(2*u)+Crs*sin(2*u))*niuDOT2;
% iDOT2=IDOT+2*(Cic*cos(2*u)+Cis*sin(2*u))*niuDOT2;

% coordinate in orbit plane coordinates
xp = r*cos(u);
yp = r*sin(u);
rc = [xp;yp;0];

xpDOT = rDOT*cos(u) - yp*uDOT;
ypDOT = rDOT*sin(u) + xp*uDOT;
rcDOT = [xpDOT; ypDOT;0];

R_ie=angle2dcm( 0 , -i, -OMEGA, 'ZXZ');
R_ieDOT=[-sin(OMEGA)*OMEGADOT, -cos(OMEGA)*cos(i)*OMEGADOT+ sin(OMEGA)*sin(i)*iDOT,  cos(OMEGA)*sin(i)*OMEGADOT+sin(OMEGA)*cos(i)*iDOT; 
          cos(OMEGA)*OMEGADOT, -sin(OMEGA)*cos(i)*OMEGADOT- cos(OMEGA)*sin(i)*iDOT,  sin(OMEGA)*sin(i)*OMEGADOT-cos(OMEGA)*cos(i)*iDOT;
          0                  ,  cos(i)*iDOT                                       , -sin(i)*iDOT                  ];
% xkDOT = ( xpDOT-yp*cos(i)*OMEGADOT )*cos(OMEGA)- ( xp*OMEGADOT+ypDOT*cos(i)-yp*sin(i)*iDOT )*sin(OMEGA);
% ykDOT = ( xpDOT-yp*cos(i)*OMEGADOT )*sin(OMEGA) + ( xp*OMEGADOT+ypDOT*cos(i)-yp*sin(i)*iDOT )*cos(OMEGA);
% zkDOT = ypDOT*sin(i) + yp*cos(i)*iDOT;
% VX=[xkDOT;ykDOT; zkDOT];
VX = R_ieDOT*rc+R_ie*rcDOT;

end

