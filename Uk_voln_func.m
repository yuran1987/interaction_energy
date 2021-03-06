function y=Uk_voln_func(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a1,a2,sigma)
  F1=@(k) 16*A1./(4+a1.^2.*k.^2).^2.+B1*(1.0-3*a1.^2.*k.^2+2*a1.^4.*k.^4)./(a1.^2.*k.^2+1.0).^4 .+C1*(1.0-5*a1.^2.*k.^2)./(a1.^2.*k.^2+1.0).^4.+D1./(1.0+a1.^2.*k.^2).^3;
  F2=@(k) 16*A2./(4+a2.^2.*k.^2).^2.+B2.*(1.0-3*a2.^2.*k.^2+2*a2.^4.*k.^4)./(a2.^2.*k.^2+1.0).^4 .+C2.*(1.0-5*a2.^2.*k.^2)./(a2.^2.*k.^2+1.0).^4.+D2./(1.0+a2.^2.*k.^2).^3;
  y = (1.0-F1(k)./Z1).^2.*(1-F2(k)./Z2).^2.*exp(-k.^2.*sigma.^2)./k;
endfunction
