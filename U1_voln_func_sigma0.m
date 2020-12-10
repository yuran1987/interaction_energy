function y=U1_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a)

e2=14.4;
R1=-(A1 + B1 + C1 + D1 - Z1)/Z1;
R2=(exp(-r/a)*(8*a^3*(B1+C1+D1)+a^2*(6*B1+8*C1+5*D1)*r+a*(2*B1+4*C1+D1)*r^2+(B1+C1)*r^3))/(16*a^3*Z1);

R3=(A1*exp(-2*r/a)*(a+r))/(2*a*Z1);

if r>0
  y=Z1*e2*(R1+2*R2+2*R3)/r;
else
  y=0;
endif
endfunction






















