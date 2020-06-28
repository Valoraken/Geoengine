%% Task 3 
% a) elevation = 60
A_a=[-sind(60) 0 -cosd(60) 1;-sind(60) -cosd(60) 0 1;-sind(60) 0 cosd(60) 1;-sind(60) cosd(60) 0 1;-1 0 0 1]
Qxx_a=inv(A_a'*A_a)
GDOP_a=sqrt(trace(Qxx_a))
PDOP_a=sqrt(Qxx_a(1,1)+Qxx_a(2,2)+Qxx_a(3,3));
HDOP_a=sqrt(Qxx_a(2,2)+Qxx_a(3,3));% local: Z->n Y->e X->u
VDOP_a=sqrt(Qxx_a(1,1));
TDOP_a=sqrt(Qxx_a(4,4));
DOP_a=[GDOP_a PDOP_a HDOP_a VDOP_a TDOP_a]
Co_a=Qxx_a(1,4)/sqrt(Qxx_a(1,1)*Qxx_a(4,4))
% b)
A_b=[-sind(45) 0 -cosd(45) 1;-sind(45) -cosd(45) 0 1;-sind(45) 0 cosd(45) 1;-sind(45) cosd(45) 0 1;-1 0 0 1]
Qxx_b=inv(A_b'*A_b)
GDOP_b=sqrt(trace(Qxx_b));
PDOP_b=sqrt(Qxx_b(1,1)+Qxx_b(2,2)+Qxx_b(3,3));
HDOP_b=sqrt(Qxx_b(2,2)+Qxx_b(3,3));% local: Z->n Y->e X->u
VDOP_b=sqrt(Qxx_b(1,1));
TDOP_b=sqrt(Qxx_b(4,4));
DOP_b=[GDOP_b PDOP_b HDOP_b VDOP_b TDOP_b]
Co_b=Qxx_b(1,4)/sqrt(Qxx_b(1,1)*Qxx_b(4,4))
% c)
A_c=[-sind(30) 0 -cosd(30) 1;-sind(30) -cosd(30) 0 1;-sind(30) 0 cosd(30) 1;-sind(30) cosd(30) 0 1;-1 0 0 1]
Qxx_c=inv(A_c'*A_c)
GDOP_c=sqrt(trace(Qxx_c));
PDOP_c=sqrt(Qxx_c(1,1)+Qxx_c(2,2)+Qxx_c(3,3));
HDOP_c=sqrt(Qxx_c(2,2)+Qxx_c(3,3));% local: Z->n Y->e X->u
VDOP_c=sqrt(Qxx_c(1,1));
TDOP_c=sqrt(Qxx_c(4,4));
DOP_c=[GDOP_c PDOP_c HDOP_c VDOP_c TDOP_c]
Co_c=Qxx_c(1,4)/sqrt(Qxx_c(1,1)*Qxx_c(4,4))
% d)
A_d=[-sind(0) 0 -cosd(0) 1;-sind(0) -cosd(0) 0 1;-sind(0) 0 cosd(0) 1;-sind(0) cosd(0) 0 1;-1 0 0 1]
Qxx_d=inv(A_d'*A_d)
GDOP_d=sqrt(trace(Qxx_d));
PDOP_d=sqrt(Qxx_d(1,1)+Qxx_d(2,2)+Qxx_d(3,3));
HDOP_d=sqrt(Qxx_d(2,2)+Qxx_d(3,3));% local: Z->n Y->e X->u
VDOP_d=sqrt(Qxx_d(1,1));
TDOP_d=sqrt(Qxx_d(4,4));
DOP_d=[GDOP_d PDOP_d HDOP_d VDOP_d TDOP_d]
Co_d=Qxx_d(1,4)/sqrt(Qxx_d(1,1)*Qxx_d(4,4))
% e)
A_e=[-sind(0) -cosd(45) -cosd(0)*cosd(45) 1;-sind(0) -cosd(0)*cosd(45) cosd(45) 1;-sind(0) cosd(45) cosd(0)*cosd(45) 1;-sind(0) cosd(0)*cosd(45) -cosd(45) 1;-1 0 0 1]
Qxx_e=inv(A_e'*A_e)
GDOP_e=sqrt(trace(Qxx_e));
PDOP_e=sqrt(Qxx_e(1,1)+Qxx_e(2,2)+Qxx_e(3,3));
HDOP_e=sqrt(Qxx_e(2,2)+Qxx_e(3,3));% local: Z->n Y->e X->u
VDOP_e=sqrt(Qxx_e(1,1));
TDOP_e=sqrt(Qxx_e(4,4));
DOP_e=[GDOP_e PDOP_e HDOP_e VDOP_e TDOP_e]
Co_e=Qxx_e(1,4)/sqrt(Qxx_e(1,1)*Qxx_e(4,4))