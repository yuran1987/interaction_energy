function y=F_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,sigma)%первая производная
if sigma>0
  e2=14.4;
  k_max = inf;
  dU=(4*pi)^2*Z1*Z2*e2*integral(@(k) Uk(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a1,a2,sigma).*(k.*cos(k.*r) - sin(k.*r)./r),0,k_max,'RelTol',0,'AbsTol',1e-12)./((2*pi)^3*r);
else
  f = @(r) U4_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,0);
  h = 0.001;
  dU = (f(r+h)-f(r-h))/(2*h);
endif
  y = -dU;
endfunction