function y=U1_molier(Z1,Z2,r,a)
e2=14.4;
R1=0.05*exp(-6*r/a);
R2=0.175*exp(-0.3*r/a);
R3=0.275*exp(-1.2*r/a);

if r>0
  y=2*Z1*Z2*e2*(R1+R2+R3)/(r);
else
  y=0;
endif
endfunction