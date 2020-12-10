function y=U2_molier(Z1,Z2,r,a)
e2=14.4;
R1=(exp(-0.3*r/a)*(0.0483289*a-0.0091875*r))/a;
R2=(exp(-1.2*r/a)*(0.354292*a-0.09075*r))/a;
R3=(exp(-6*r/a)*(0.0973794*a-0.015*r))/a;

if r>0
  y=2*Z1*Z2*e2*(R1+R2+R3)/(r);
else
  y=0;
endif
endfunction
