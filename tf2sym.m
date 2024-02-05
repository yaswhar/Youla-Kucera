function Gsym = tf2sym(Gtf)
% Transform transfer functions to symbolic expressions
[num, den] = tfdata(Gtf, 'v');
syms s
symNum = poly2sym(num,s);
symDen = poly2sym(den,s);
Gsym = symNum / symDen;
end