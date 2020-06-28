% Lab 4 part 2
% Xiao, Tianqi 3371477
% 08.01.19

function coe_PLx = polyLagrange(x,y)
% This function is to compute the coefficients of the Lagrange
% interpolating polynomial PL(x) in ascending power of x for a given data points (xi,yi)
%   INPUT:
% x is the x coordinates of the data points
% y is the y coordinates of the data points
%   OUTPUT:
% coe_PLx is the coefficients of the Lagrange interpolating polynomial PL(x)
% the length of x and y must be the same
if length(x)~=length(y)
    error('the length of the inputs should be the same')
end

n=length(x)
for i=1:n
    lix=1;
    for j = 1:n
        if j~=i
            % descending order of x
            lix=conv([-x(j)/(x(i)-x(j)),1/(x(i)-x(j))],lix);
        else
            lix=lix*1;
        end
    end
    Lx(i,:)=lix;
end
coe_PLx=y'*Lx;
end

