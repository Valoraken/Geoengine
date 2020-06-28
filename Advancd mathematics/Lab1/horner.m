% Xiao Tianqi 3371477

function OutputValue = horner(an,x0)
% This funtion is for caculating the value of the polynomial
%      Pn(x)=an*x^n+a(n-1)*x^(n-1)+...+a1x+a0 with degree n by horner scheme
% Input: an is the coefficient vector
%        x0 is a guess root
%
% Output: OutputValue is the value of Pn(x)

OutputValue=0;
n=length(an);
% if i=n,bi=an
% else bi=ai+b(i+1)*x0
b=zeros(1,n);
b(n)= an(n);
for i=n-1:-1:1
   
   b(i)=an(i)+b(i+1)*x0;

end
%present the Horner scheme
  fprintf('x0=%d\n',x0)
  fprintf('an:\n')
  fprintf('%d\t\t',an);
  fprintf('\nbn:\n')
  fprintf('\t%d*x0 ',b(2:n));
  fprintf('\n-----------------------------------\n');
  fprintf('b0=P(x0)=%d\n\n',b(1));

  OutputValue=b(1);
    
end

