y1=[2;4;5];
A1=[9 2 10;4 -6 1;3 1 77];
rank(A1)
Ay1=[A1 y1];
rank(Ay1)
x1=A1\y1
y2=[5.2 0 1.1]';
A2=[-6 2;8 6;1 3];
rank(A2)
Ay2=[A2 y2];
rank(Ay2)
y3=[-18.4 -16.2 -12.6]';
A3=[-6 2;8 6;1 3];
rank(A3)
Ay3=[A3 y3];
rank(Ay3)
x3=A3\y3
y4=[5 -10 1]';
A4=[-2.5 1.2 0.5;5 -2.4 -1 ;-0.5 -2.7 77];
rank(A4)
Ay4=[A4 y4];
rank(Ay4)
x4=A4\y4
y5=[1.4 -2.7 4]';
A5=[-2.5 1.2 0.5;5 -2.4 -1 ;-0.5 -2.7 3.3];
rank(A5)
Ay5=[A5 y5];
rank(Ay5)
x5=A5\y5
y6=[8 3]';
A6=[1 2 3 4;5 6 7 8];
rank(A6)
Ay6=[A6 y6];
rank(Ay6)
x6=A6\y6
y7=[-2 4 5]';
A7=[9 2 10;4 -6 1;3 1 3.5];
rank(A7)
Ay7=[A7 y7];
rank(Ay7)
x7=A7\y7
y8=[-1.4 -0.1 5.4]';
A8=[6.25 -1 77;-1 0.16 -1.28;8 -1.28 10.24];
rank(A8)
Ay8=[A8 y8];
rank(Ay8)
x8=A8\y8