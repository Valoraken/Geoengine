clear

%%  Task 2
syms t
t_p=90;
t_q=0:90;
psi_pq=sin(deg2rad(t_p))*sin(deg2rad(t_q))+cos(deg2rad(t_p))*cos(deg2rad(t_q)); 

for l=0:100
    
    pl_psipq=1/(2^l*factorial(l))*diff((t^2-1)^l,t,l);
  
    P_left(l+1,:)=double(subs(pl_psipq,t,psi_pq));

for m=0:l
P_lm=zeros(size(psi_pq));
if m==0
    N=sqrt(2*l+1);
else
    N=sqrt(2*(2*l+1)*(factorial(l-m)/factorial(l+m)));
end
P_lm=P_lm+double(subs((1-t^2).^(m/2).*(diff(pl_psipq,t))*N,t,cos(deg2rad(t_p)))).*double(subs((1-t^2).^(m/2).*(diff(pl_psipq,t))*N,t,cos(deg2rad(t_q))));

end

P_right(l+1,:)=1/(2*l+1).*P_lm;
l
end