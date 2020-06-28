syms a b
R1x=[1 0 0; 0 cos(a) sin(a);0 -sin(a) cos(a)]
R1y=[cos(a) 0 sin(a); 0 1 0;-sin(a) 0 cos(a)]
R1z=[cos(a) sin(a) 0;-sin(a) cos(a) 0 ;0 0 1]
R2x=[1 0 0; 0 cos(b) sin(b);0 -sin(b) cos(b)]
R2y=[cos(b) 0 sin(b); 0 1 0;-sin(b) 0 cos(b)]
R2z=[cos(b) sin(b) 0;-sin(b) cos(b) 0 ;0 0 1]

R1x*R1y*R1z*R2x*R2y*R2z-R1x*R2x*R1y*R2y*R1z*R2z