% Xiao, Tianqi 3371477
% Lab 5 27.11.18

function [A] = jacobipq(A,iterMax,tol)
% This function is to find out symmetric eigen values iteratively by (p,q)rotation of Jacobi
%                         A[L+1]=U'*A[L]*U
% with U(i,i)=1;  U(p,p)=U(q,q)=cos(phi);   U(p,q)=sin(phi); U(q,p)=-sin(phi)
% with cot(2*phi) = (A(q,q)-A(p,p))/(2*A(p,q))
%
% IMPUT:
% A is a symmetric matrix
% iterMax is the maximum interations
% tol is the tolerance, the values smaller than tol will be set to zero
%
% OUTPUT:
% A is the diagonal matrix of the input matrix or the result after the maximun iterations

%% check input
if A'~=A
    error('A is not symmetric!');
end
%% initialize
count = 1;
U=eye(size(A));
Amax=0;
p=1; q=1;
%% interation
while count < iterMax
    %find out the off-diagonal element A(p,q) with largest absolute value
    for r=1:length(A)
        for c=1:length(A)
            if r~=c && Amax<abs(A(r,c))
                Amax=abs(A(r,c));
                p=r;
                q=c;
            end
        end
    end
    % if A is diagonal, stop the interation
    if Amax==0
        break;
    end
    phi=1/2*acot((A(q,q)-A(p,p))/(2*A(p,q)));
    % build U
    U(p,p)=cos(phi);
    U(q,q)=cos(phi);
    U(p,q)=sin(phi);
    U(q,p)=-sin(phi);
    
    A=U'*A*U;
    % set small values to zero
    A(A<tol)=0;
    % show the elements of the first step
    if count==1
        A
    end
    count=count+1;

    % reset U and Amax
    U=eye(size(A));
    Amax = 0;
end
end

