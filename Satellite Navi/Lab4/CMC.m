function CMC = CMC(f1,f2,phi1,phi2,r1,r2)
   c=physconst('Lightspeed')
PhiIF= (f1.^2.*phi1-f2.^2.*phi2)/(f1.^2-f2.^2);
rIF=(f1.^2.*r1-f2.^2.*r2)/(f1.^2-f2.^2);
n1=f1.^2/(f1.^2-f2.^2);
n2=-f2.^2/(f1.^2-f2.^2);
LambdaIF=c/(n1*f1+n2*f2);
CMC=rIF-PhiIF*LambdaIF;
end