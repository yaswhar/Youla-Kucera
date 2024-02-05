function Gtf = sym2tf(Gsym)
% Transform symbolic expressions to transfer functions
[symNum, symDen] = numden(Gsym);
TFnum = sym2poly(symNum);  
TFden = sym2poly(symDen);
Gtf = tf(TFnum, TFden);
end