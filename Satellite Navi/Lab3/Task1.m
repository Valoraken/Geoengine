% Satellite Navigation
% Tianqi Xiao
% 3371477

% Task 1
%% GPS C/A code
ini_GPS= [1 1 1 1 1 1 1 1 1 1];
% Generate G1 code
shifted=[0 0 0 0 0 0 0 0 0 0];
temp=ini_GPS;
output_G1=[];
while isequal(shifted,ini_GPS)==0
    shifted=[xor(temp(3), temp(10)) temp(1:end-1)];
    temp= shifted;
    output_G1=[output_G1;shifted];
end
% pick 10th number as G1 code 
G1=output_G1(:,10)';
% Generate G2 code
shifted=[0 0 0 0 0 0 0 0 0 0];
temp=ini_GPS;
output_G2=[];
while isequal(shifted,ini_GPS)==0
    shifted=[xor(temp(10),xor(temp(9),xor(temp(8),xor(temp(6),xor(temp(2), temp(3)))))) temp(1:end-1)];
    temp= shifted;
    output_G2=[output_G2;shifted];
end
% pick xor(i,j) number as G2 code 
% C/A=xor(G1,G2)
%% GLONASS PRN
% Generate shift register
% G(X)=1+X(5)+X(9)
ini_GLONASS=[1 1 1 1 1 1 1 1 1];
shifted=[0 0 0 0 0 0 0 0 0];
temp=ini_GLONASS;
output_GLONASS=[];
while isequal(shifted,ini_GLONASS)==0
    shifted=[xor(temp(5), temp(9)) temp(1:end-1)];
    temp= shifted;
    output_GLONASS=[output_GLONASS;shifted];
end
% pick 7th number as PRN code
PRN=output_GLONASS(:,7);
PRN_20=[PRN(1:20)';PRN(end-19:end)']
%%  ±±¶·Beidou B1I ranging code
ini_beidou=[0 1 0 1 0 1 0 1 0 1 0];

% Generate G1 code
shifted=[0 0 0 0 0 0 0 0 0 0 0];
temp=ini_beidou;
output_BG1=[];
while isequal(shifted,ini_beidou)==0
    shifted=[xor(temp(11),xor(temp(10),xor(temp(9),xor(temp(8),xor(temp(1), temp(7)))))) temp(1:end-1)];
    temp= shifted;
    output_BG1=[output_BG1;shifted];
end
BG1=output_BG1(1:end-1,11);
% Generate G2 code
shifted=[0 0 0 0 0 0 0 0 0 0 0];
temp=ini_beidou;
output_BG2=[];
while isequal(shifted,ini_beidou)==0
    shifted=[xor(temp(11),xor(temp(9),xor(temp(8),xor(temp(5),xor(temp(4),xor(temp(3),xor(temp(1), temp(2)))))))) temp(1:end-1)];
    temp= shifted;
    output_BG2=[output_BG2;shifted];
end

% sequence of G2 code for satellit1 1-15
seq_G2=[1 3;1 4;1 5;1 6;1 8;1 9;1 10;1 11;2 7;3 4;3 5;3 6;3 8;3 9;3 10];
BG2=xor(output_BG2(1:end-1,seq_G2(:,1)),output_BG2(1:end-1,seq_G2(:,2)));
B1I=1*xor(BG2,BG1)';
B1I_10=[B1I(:,1:10);B1I(:,end-9:end)]


data=load('bds.txt');
figure 
hold on
index=[];
for i=1:15
code=B1I(i,:);
    CC=zeros(1,length(data));
    for j=0:length(data)-1
        xi=circshift(data,j);
        CC(j+1)=sum(xi.*code);
    end
    plot(0:length(data)-1,CC)
if max(abs(CC))>100
    index=[index; i find(abs(CC)>100) CC(find(abs(CC)>100))];
end
end

legend(num2str((1:15)'));
