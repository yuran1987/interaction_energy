function y=ddU_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,sigma)%вторая производная
if sigma>0
  e2=14.4;
  k_max = inf;
  ddsin = @(k) 2.*sin(k.*r)./r.^2 - k.^2.*sin(k.*r) - 2*k.*cos(k.*r)./r;
  y=(4*pi)^2*Z1*Z2*e2*integral(@(k) Uk(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a1,a2,sigma).*ddsin(k),0,k_max,'RelTol',0,'AbsTol',1e-12)./((2*pi)^3*r);
else
  f = @(r) U4_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,0);
  h = 0.001;
  y = (f(r+h)-2*f(r)+f(r-h))/(h^2);
endif
endfunction