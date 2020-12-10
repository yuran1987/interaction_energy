function y=U4_voln_func(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1,a2,sigma, is_hartree_fock)
  e2=14.4;
  k_max = inf;
if is_hartree_fock == 1 %через волновые функции Хартри-Фока
    y = (4*pi)^2*Z1*Z2*e2*integral(@(k) Uk_voln_func_over_hartreefock(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,sigma).*sin(k.*r),0,k_max,'RelTol',0,'AbsTol',1e-12)./((2*pi)^3*r);
else %через решение уравнения Шрёдингера  
    if sigma>0
      y=(4*pi)^2*Z1*Z2*e2*integral(@(k) Uk_voln_func(k,A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,a1,a2,sigma).*sin(k.*r),0,k_max,'RelTol',0,'AbsTol',1e-12)./((2*pi)^3*r);
    else
      y = U4_voln_func_sigma0(A1,B1,C1,D1,A2,B2,C2,D2,Z1,Z2,r,a1);
    endif
endif  
endfunction






















