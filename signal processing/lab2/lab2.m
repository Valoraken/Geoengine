clear
M1=[1 1 1 9 9 9 1 1 1 9 9 9 1 1 1;
    1 1 1 9 9 9 1 1 1 9 9 9 1 1 1;
    1 1 1 9 9 9 1 1 1 9 9 9 1 1 1;
    9 9 9 9 9 9 9 9 9 9 9 9 9 9 9;
    9 9 9 9 9 9 9 9 9 9 9 9 9 9 9;
    9 9 9 9 9 9 9 9 9 9 9 9 9 9 9;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]
Sizem1=size(M1);
M2=M1
M4=M1
M3=28.*M1;
figure(1)

imshow(M3)
%model1 
%a=[1 1 1 ;2 2 2; 1 2 1]
a=[-1 0 1;-2 0 2;-1 0 1];
b=[1 2 1; 0 0 0;-1 -2 -1];
for i=2:Sizem1(1)-1
for j=2:Sizem1(2)-1
A=double(M1(i-1:i+1,j-1:j+1)).*a;
M2(i,j)=sum(sum(A));
figure(2)
imshow(M2)
end
end

for i=2:Sizem1(1)-1
for j=2:Sizem1(2)-1
B=double(M1(i-1:i+1,j-1:j+1)).*b;
M4(i,j)=sum(sum(B));
figure(3)
imshow(M4)
end
end
M2
M4